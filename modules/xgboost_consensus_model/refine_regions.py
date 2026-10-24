"""
refine_regions.py -- turn scored consensus candidates into the FINAL region set.

Consensus candidates come out of _merge_candidates() with UNION coordinates, so
a single segmental call (e.g. a multi-megabase control-FREEC or cnvkit segment)
can bridge exon-level calls that are hundreds of kb apart into one region. Two
things then go wrong downstream if that region is used as-is:

  1. The region is annotated as a whole, so its gene list / classification /
     score describe the union span rather than the places evidence exists.
  2. Fanning that span out across the capture BED invents one row per panel
     exon it touches, including exons where no caller called anything.

This script therefore does the breakdown HERE, before annotation, and keeps only
exons that actually carry caller support:

  * intersect each candidate with the capture BED,
  * drop exons not overlapped by any non-outlier member call,
  * re-merge surviving exons of the same sample/chr/type across gaps <= --gap.

Running before CLASSIFY_AND_ANNOTATE is the point: the region set is final by
the time anything annotates it, so genes, Classification, Total score and the
dosage-sensitive gene list are all computed on the coordinates they describe.

"Outlier" reuses MAX_CLUSTER_SIZE_RATIO from cnv_consensus_model, with the same
meaning it has in _merge_candidates: a member call more than that many times the
median member-call size contributes support nowhere, because it is the thing
doing the spurious bridging.

Input and output are both the score_samples_xgb.py column set:
    sample | chr | start | end | type | n_callers | callers_supporting |
    consensus_score | cn_label
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from cnv_consensus_model import (ALL_CALLERS, MAX_CLUSTER_SIZE_RATIO,  # noqa: E402
                                 load_sample_calls, norm_chr)

# A member call spanning exons of more than this many distinct panel genes is
# treated as a segmental/blanket call: it supports the candidate's existence but
# cannot say which exons are affected. 2 allows a genuine whole-gene deletion,
# and a call straddling two adjacent panel genes, while rejecting the
# multi-megabase control-FREEC/cnvkit segments that do the spurious bridging.
MAX_SUPPORT_GENES = 2

COLUMNS = ["sample", "chr", "start", "end", "type", "n_callers",
           "callers_supporting", "consensus_score", "cn_label"]


def load_panel(bed_path):
    """{chr: [(start, end, name)]} from a BED, converted to 1-based inclusive
    starts to match the candidate coordinates. `name` is BED column 4 when
    present (the capture BEDs here carry the gene symbol) else None."""
    panel = {}
    with open(bed_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            try:
                start, end = int(f[1]), int(f[2])
            except ValueError:
                continue
            name = f[3].strip() if len(f) >= 4 and f[3].strip() else None
            panel.setdefault(norm_chr(f[0]), []).append((start + 1, end, name))
    for c in panel:
        panel[c].sort()
    return panel


def _genes_touched(interval, exons):
    """Distinct panel gene names whose exons the interval overlaps."""
    s, e = interval
    return {n for xs, xe, n in exons if n and xs < e and xe > s}


def support_intervals(calls_by_caller, row, exons):
    """Member calls that may validate an exon, with the bridging calls removed.

    A member call is one from a supporting caller, on the same chr and of the
    same type, overlapping the candidate. Two filters then apply, because a
    call that blankets the whole neighbourhood is evidence that *something*
    is there, not evidence about any particular exon:

      * scale -- more than MAX_CLUSTER_SIZE_RATIO times the SMALLEST member
        call. Using the smallest rather than the median matters: with only two
        members `sizes[len//2]` is the larger of the pair, which makes the test
        vacuous exactly when a small exon call is paired with a big segment.
      * breadth -- spanning exons of more than MAX_SUPPORT_GENES distinct panel
        genes. This is what catches a lone multi-megabase segmental call, which
        no relative-scale test can flag because it is its own scale reference.
        Skipped when the BED carries no gene names.

    Returns [] when nothing survives; the caller then leaves the region intact
    rather than fanning it across exons nothing actually called.
    """
    supporting = set(str(row["callers_supporting"]).split(","))
    members = []
    for caller, calls in calls_by_caller.items():
        if caller not in supporting:
            continue
        for c in calls:
            if (norm_chr(c["chr"]) == row["chr"] and c["type"] == row["type"]
                    and c["start"] < row["end"] and c["end"] > row["start"]):
                members.append((c["start"], c["end"]))
    if not members:
        return []

    smallest = min(e - s for s, e in members)
    threshold = smallest * MAX_CLUSTER_SIZE_RATIO
    kept = [m for m in members if (m[1] - m[0]) <= threshold]
    kept = [m for m in kept
            if len(_genes_touched(m, exons)) <= MAX_SUPPORT_GENES
            or not any(n for _, _, n in exons)]
    return kept


def refine(scored, calls_by_caller, panel):
    """One row per supported panel exon. Candidates with no panel overlap, or
    no surviving support, are passed through unchanged rather than dropped --
    losing a called region silently is worse than reporting a coarse one."""
    out = []
    for _, row in scored.iterrows():
        chrom_exons = panel.get(row["chr"], [])
        exons = [(s, e, n) for s, e, n in chrom_exons
                 if s < row["end"] and e > row["start"]]
        support = support_intervals(calls_by_caller, row, chrom_exons)
        if not exons or not support:
            out.append({**dict(row), "_broken_down": False})
            continue
        kept = [(s, e) for s, e, _ in exons
                if any(cs < e and ce > s for cs, ce in support)]
        if not kept:
            out.append({**dict(row), "_broken_down": False})
            continue
        for s, e in kept:
            out.append({**dict(row), "start": s, "end": e, "_broken_down": True})
    return pd.DataFrame(out, columns=COLUMNS + ["_broken_down"])


def smooth(df, gap):
    """Merge same sample/chr/type rows separated by <= gap bp.

    Only rows that came out of the exon breakdown take part. A region passed
    through unrefined keeps a caller-scale span, which can be many megabases,
    and letting one into the running `cur["end"]` would make it absorb every
    later row on the chromosome -- the same snowballing that _merge_candidates
    does with a bridging segmental call, reintroduced one stage later. Those
    rows are emitted exactly as scored instead.
    """
    if df.empty:
        return df.drop(columns=["_broken_down"], errors="ignore")

    passthrough = df[~df["_broken_down"]].drop(columns=["_broken_down"])
    breakable = (df[df["_broken_down"]].drop(columns=["_broken_down"])
                 .sort_values(["sample", "chr", "type", "start"])
                 .reset_index(drop=True))

    segments, cur = [], None
    for _, r in breakable.iterrows():
        if (cur is not None and r["sample"] == cur["sample"]
                and r["chr"] == cur["chr"] and r["type"] == cur["type"]
                and r["start"] - cur["end"] <= gap):
            cur["end"] = max(cur["end"], r["end"])
            cur["n_callers"] = max(cur["n_callers"], r["n_callers"])
            cur["consensus_score"] = max(cur["consensus_score"], r["consensus_score"])
            cur["callers_supporting"] = ",".join(sorted(
                set(str(cur["callers_supporting"]).split(","))
                | set(str(r["callers_supporting"]).split(","))))
            cur["cn_label"] = ";".join(dict.fromkeys(
                str(cur["cn_label"]).split(";") + str(r["cn_label"]).split(";")))
        else:
            if cur is not None:
                segments.append(cur)
            cur = dict(r)
    if cur is not None:
        segments.append(cur)

    merged = pd.DataFrame(segments, columns=COLUMNS)
    return (pd.concat([merged, passthrough], ignore_index=True)
            .sort_values(["sample", "chr", "start", "type"])
            .reset_index(drop=True)[COLUMNS])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--scored", required=True, help="Scored candidate TSV")
    ap.add_argument("--sample_dir", required=True, help="This sample's caller files")
    ap.add_argument("--panel_bed", required=True, help="Capture/panel BED")
    ap.add_argument("--gap", type=int, default=10000, help="Smoothing gap in bp")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    scored = pd.read_csv(args.scored, sep="\t", dtype={"chr": str, "sample": str})
    if scored.empty:
        scored.to_csv(args.out, sep="\t", index=False)
        print("refine_regions: no scored candidates, nothing to refine")
        return

    scored["chr"] = scored["chr"].map(norm_chr)
    panel = load_panel(args.panel_bed)
    calls_by_caller = load_sample_calls(Path(args.sample_dir), ALL_CALLERS)

    refined = smooth(refine(scored, calls_by_caller, panel), args.gap)
    refined.to_csv(args.out, sep="\t", index=False)

    print(f"refine_regions: {len(scored)} candidate region(s) "
          f"(total span {int((scored['end'] - scored['start']).sum()):,} bp) "
          f"-> {len(refined)} final region(s) "
          f"(total span {int((refined['end'] - refined['start']).sum()):,} bp)")


if __name__ == "__main__":
    main()
