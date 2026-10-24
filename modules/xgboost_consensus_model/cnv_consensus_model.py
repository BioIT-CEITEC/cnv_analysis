"""
CNV caller parsing, overlap logic and candidate-region merging.
==============================================================
Shared library for the XGBoost consensus model. It owns everything that turns
raw per-caller output into scored-candidate regions EXCEPT the classifier
itself:

  * one parser per caller, normalising to {chr, start, end, type, score, cn}
  * overlap/coverage helpers (MIN_OVERLAP_FRAC)
  * _merge_candidates(): clusters overlapping calls of the same chr+type into
    the candidate regions that get scored, and decides their coordinates
  * ground-truth loading, used by the offline training/evaluation scripts

The caller parsers, ground-truth loader and overlap logic are inherited
verbatim from the original reference model, so behaviour stays comparable
across model generations.

Consumers: score_samples_xgb.py (this module), modelling/merge_from_bronco.py,
modelling/src/{events,eda}.py.
"""

import re
import warnings
from pathlib import Path

import pandas as pd

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Every caller the parsers understand.
ALL_CALLERS = ["cnvkit", "freec", "exomeDepth", "panelcnMOPS", "cnMOPS",
               "XHMM", "conifer", "gatk"]

# Minimum fraction of the SMALLER interval that must be covered for a match.
# (Kept identical to the reference: a large caller segment containing a small
#  GT exon still counts as a hit.)
MIN_OVERLAP_FRAC = 0.5

# When building merged candidates, any call whose size exceeds
# MAX_CLUSTER_SIZE_RATIO × the median call size in the cluster is an outlier
# (e.g. conifer reporting a gene-level deletion while exon-level callers report
# a single exon).  Outlier calls are counted in n_callers but their coordinates
# are excluded from the candidate boundaries so they cannot inflate a small
# exon-level candidate to a multi-megabase region.
MAX_CLUSTER_SIZE_RATIO = 50


# ---------------------------------------------------------------------------
# Utility helpers  (verbatim from reference)
# ---------------------------------------------------------------------------

def norm_chr(c):
    """Normalise chromosome name to 'chrN' form."""
    c = str(c).strip()
    if not c.startswith("chr"):
        c = "chr" + c
    return c


def norm_type(t):
    """Return 'DEL' or 'DUP' or None."""
    if pd.isna(t):
        return None
    t = str(t).lower()
    if any(k in t for k in ["del", "loss", "cn0", "cn1"]):
        return "DEL"
    if any(k in t for k in ["dup", "gain", "amp"]):
        return "DUP"
    m = re.match(r"cn(\d+)", t)
    if m:
        cn = int(m.group(1))
        if cn < 2:
            return "DEL"
        if cn > 2:
            return "DUP"
    return None


def fmt_cn(v):
    """Render one caller's copy number for cn_label.

    Two pandas artefacts have to be absorbed here. The "cn" column is built
    from call.get("cn"), and cnvkit/freec/cnMOPS never set that key while
    exomeDepth/conifer set it to None -- so the column contains None values,
    which pandas converts to float NaN when the candidate frame is built.
    NaN is not None, so an `is not None` test lets it through and it renders
    as the literal "nan". That same NaN also promotes the column to float64,
    so an integer CN of 3 arrives here as 3.0.
    """
    if v is None or pd.isna(v):
        return "NA"
    if isinstance(v, float) and v.is_integer():
        return str(int(v))
    return str(v)


def gt_coverage(s1, e1, s2, e2):
    """Fraction of the SMALLER interval covered by its overlap with the larger."""
    overlap = max(0, min(e1, e2) - max(s1, s2))
    if overlap == 0:
        return 0.0
    min_len = min(e1 - s1, e2 - s2)
    if min_len <= 0:
        return 0.0
    return overlap / min_len


def overlaps(s1, e1, s2, e2, min_frac=0.0):
    """True if the two intervals overlap with min-length coverage >= min_frac."""
    return gt_coverage(s1, e1, s2, e2) >= min_frac


# ---------------------------------------------------------------------------
# Ground-truth loading  (Format B: *_selected_regions.tsv)
# ---------------------------------------------------------------------------

def merge_gt_variants(variants):
    """
    Deduplicate and merge exon-level GT entries into variant-level events.
    Two entries merge when on the same chr, same type, and start within 10 kb
    of the current end (consecutive panel exons of the same gene).
    """
    if not variants:
        return []
    GAP = 10_000
    seen, unique = set(), []
    for v in variants:
        key = (v["chr"], v["start"], v["end"], v["type"])
        if key not in seen:
            seen.add(key)
            unique.append(v)
    unique.sort(key=lambda x: (x["chr"], x["type"], x["start"]))
    merged = []
    cur = dict(unique[0])
    for v in unique[1:]:
        if (v["chr"] == cur["chr"] and v["type"] == cur["type"]
                and v["start"] <= cur["end"] + GAP):
            cur["end"] = max(cur["end"], v["end"])
        else:
            merged.append(cur)
            cur = dict(v)
    merged.append(cur)
    return merged


def load_ground_truth_dir(gt_dir):
    """
    Directory of per-sample TSVs named <SAMPLE>_selected_regions.tsv.
    Columns: sample, chrom, start, end, name, variant_type, type
    'type' holds DEL / DUP directly.
    Returns (gt: dict sample -> [ {chr,start,end,type} ], all_samples: set).
    """
    gt, all_samples = {}, set()
    for tsv in sorted(Path(gt_dir).glob("*_selected_regions.tsv")):
        if tsv.name.startswith(".") or tsv.stat().st_size == 0:
            continue
        df = pd.read_csv(tsv, sep="\t", encoding="utf-8", encoding_errors="replace")
        df.columns = df.columns.str.strip()
        for _, row in df.iterrows():
            sample = str(row["sample"]).strip()
            all_samples.add(sample)
            vtype = str(row["type"]).strip().upper()
            if vtype not in ("DEL", "DUP"):
                continue
            gt.setdefault(sample, []).append({
                "chr":   norm_chr(row["chrom"]),
                "start": int(row["start"]),
                "end":   int(row["end"]),
                "type":  vtype,
            })
    gt = {s: merge_gt_variants(v) for s, v in gt.items()}
    return gt, all_samples


# ---------------------------------------------------------------------------
# Per-caller parsers  (verbatim from reference)
# Each returns a list of {chr, start, end, type, score}.
# ---------------------------------------------------------------------------

def parse_cnvkit(path):
    calls = []
    try:
        df = pd.read_csv(path, sep="\t")
        for _, r in df.iterrows():
            vtype = norm_type(r.get("SVTYPE", ""))
            if vtype is None:
                continue
            score = abs(float(r.get("FOLD_CHANGE_LOG", 0) or 0))
            calls.append({"chr": norm_chr(r["CHROM"]), "start": int(r["POS"]),
                          "end": int(r["END"]), "type": vtype, "score": score})
    except Exception:
        pass
    return calls


def parse_freec(path):
    calls = []
    try:
        df = pd.read_csv(path, sep="\t")
        for _, r in df.iterrows():
            vtype = norm_type(r.get("type", ""))
            if vtype is None:
                cn = r.get("CN", 2)
                try:
                    cn = float(cn)
                    vtype = "DEL" if cn < 2 else ("DUP" if cn > 2 else None)
                except Exception:
                    pass
            if vtype is None:
                continue
            cn = r.get("CN", 2)
            try:
                score = abs(float(cn) - 2.0)
            except Exception:
                score = 1.0
            calls.append({"chr": norm_chr(r["chromosome"]), "start": int(r["start"]),
                          "end": int(r["end"]), "type": vtype, "score": score})
    except Exception:
        pass
    return calls


def parse_exomeDepth(path):
    calls = []
    try:
        df = pd.read_csv(path, sep="\t")
        for _, r in df.iterrows():
            vtype = norm_type(r.get("type", ""))
            if vtype is None:
                continue
            bf = r.get("BF", 0)
            try:
                score = float(bf) if not pd.isna(bf) else 0.0
            except Exception:
                score = 0.0
            calls.append({"chr": norm_chr(r["chromosome"]), "start": int(r["start"]),
                          "end": int(r["end"]), "type": vtype, "score": score,
                          "cn": None})  # ExomeDepth does not report CN
    except Exception:
        pass
    return calls


def parse_panelcnMOPS(path):
    calls = []
    try:
        df = pd.read_csv(path, sep="\t")
        for _, r in df.iterrows():
            cn_raw = r.get("CN", "")
            vtype = norm_type(cn_raw)
            if vtype is None:
                continue
            try:
                rc_norm = float(r.get("RC.norm", 0) or 0)
                med_norm = float(r.get("medRC.norm", 1) or 1)
                score = abs(rc_norm / max(med_norm, 1) - 1)
            except Exception:
                score = 1.0
            # CN column is "CN0", "CN1", "CN3" etc. — extract the integer.
            cn = None
            try:
                cn_str = str(cn_raw).strip()
                if cn_str.upper().startswith("CN"):
                    cn = int(cn_str[2:])
            except Exception:
                pass
            calls.append({"chr": norm_chr(r["chromosome"]), "start": int(r["start"]),
                          "end": int(r["end"]), "type": vtype, "score": score, "cn": cn})
    except Exception:
        pass
    return calls


def parse_cnMOPS(path):
    calls = []
    try:
        df = pd.read_csv(path, sep="\t")
        if df.empty:
            return calls
        for _, r in df.iterrows():
            vtype = norm_type(r.get("CN", ""))
            if vtype is None:
                continue
            score = abs(float(r.get("mean", 0) or 0))
            calls.append({"chr": norm_chr(r.get("chromosome", r.get("seqnames"))), "start": int(r["start"]),
                          "end": int(r["end"]), "type": vtype, "score": score})
    except Exception:
        pass
    return calls


def parse_XHMM(path):
    calls = []
    try:
        df = pd.read_csv(path, sep="\t")
        for _, r in df.iterrows():
            vtype = norm_type(r.get("type", ""))
            cn_raw = r.get("CN", 2)
            if vtype is None:
                try:
                    cn_raw = float(cn_raw)
                    vtype = "DEL" if cn_raw < 2 else ("DUP" if cn_raw > 2 else None)
                except Exception:
                    pass
            if vtype is None:
                continue
            cn_raw = r.get("CN", 2)
            cn = None
            try:
                cn = int(float(cn_raw))
            except Exception:
                pass
            try:
                score = abs(float(cn_raw) - 2.0)
            except Exception:
                score = 1.0
            calls.append({"chr": norm_chr(r["chromosome"]), "start": int(r["start"]),
                          "end": int(r["end"]), "type": vtype, "score": score, "cn": cn})
    except Exception:
        pass
    return calls


def parse_conifer(path):
    calls = []
    try:
        df = pd.read_csv(path, sep="\t")
        for _, r in df.iterrows():
            vtype = norm_type(r.get("type", "")) or norm_type(r.get("CN", ""))
            if vtype is None:
                continue
            try:
                rc_norm = float(r.get("RC.norm", 0) or 0)
                med_norm = float(r.get("medRC.norm", 1) or 1)
                score = abs(rc_norm / max(med_norm, 1) - 1)
            except Exception:
                score = 1.0
            calls.append({"chr": norm_chr(r["chromosome"]), "start": int(r["start"]),
                          "end": int(r["end"]), "type": vtype, "score": score,
                          "cn": None})  # conifer CN column is always empty
    except Exception:
        pass
    return calls


def parse_gatk(path):
    """GATK gCNV VCF (gzipped/BGZF or plain). ALT <DEL>/<DUP>, END in INFO, QS score."""
    import gzip
    import io as _io
    import subprocess
    import sys
    calls = []
    path = Path(path)
    if not path.exists():
        return calls
    if str(path).endswith(".vcf.gz") and path.stat().st_size == 4096:
        sib = Path(str(path)[:-3])
        if sib.exists():
            return parse_gatk(sib)

    raw_content = None
    for opener_fn in (lambda p: gzip.open(str(p), "rb").read(),
                      lambda p: open(str(p), "rb").read()):
        try:
            data = opener_fn(path)
            if data:
                raw_content = data
                break
        except Exception:
            pass
    if raw_content is None:
        for cmd in (["zcat", str(path)], ["gunzip", "-c", str(path)]):
            try:
                result = subprocess.run(cmd, stdout=subprocess.PIPE,
                                        stderr=subprocess.DEVNULL)
                if result.stdout:
                    raw_content = result.stdout
                    break
            except FileNotFoundError:
                pass
    if raw_content is None:
        print(f"WARNING: parse_gatk could not read {path}", file=sys.stderr)
        return calls

    header = None
    for raw in _io.BytesIO(raw_content):
        line = raw.decode("utf-8", errors="replace") if isinstance(raw, bytes) else raw
        if line.startswith("##"):
            continue
        if line.startswith("#"):
            header = line.lstrip("#").strip().split("\t")
            continue
        if header is None:
            continue
        row = dict(zip(header, line.strip().split("\t")))
        vtype = norm_type(row.get("ALT", "."))
        if vtype is None:
            continue
        chrom = norm_chr(row.get("CHROM", ""))
        pos = int(row.get("POS", 0))
        end = pos
        for token in row.get("INFO", "").split(";"):
            if token.startswith("END="):
                end = int(token.split("=")[1])
                break
        fmt_keys = row.get("FORMAT", "").split(":")
        sample_col = [v for k, v in row.items()
                      if k not in ("CHROM", "POS", "ID", "REF", "ALT",
                                   "QUAL", "FILTER", "INFO", "FORMAT")]
        score = 0.0
        cn = None
        if sample_col:
            fmt = dict(zip(fmt_keys, sample_col[0].split(":")))
            try:
                score = float(fmt.get("QS", fmt.get("QA", 0)) or 0)
            except (ValueError, TypeError):
                score = 0.0
            try:
                cn = int(fmt["CN"])
            except (KeyError, ValueError, TypeError):
                cn = None
        calls.append({"chr": chrom, "start": pos, "end": end,
                      "type": vtype, "score": score, "cn": cn})
    return calls


CALLER_PARSERS = {
    "cnvkit": parse_cnvkit, "freec": parse_freec, "exomeDepth": parse_exomeDepth,
    "panelcnMOPS": parse_panelcnMOPS, "cnMOPS": parse_cnMOPS, "XHMM": parse_XHMM,
    "conifer": parse_conifer, "gatk": parse_gatk,
}

CALLER_FILE_PATTERNS = {
    "cnvkit": "*_cnvkit.tsv", "freec": "*_freec.tsv", "exomeDepth": "*_ExomeDepth.tsv",
    "panelcnMOPS": "*_panelcnMOPS.tsv", "cnMOPS": "*_cnMOPS.tsv", "XHMM": "*_xhmm.tsv",
    "conifer": "*_conifer.tsv", "gatk": "*_gatk.vcf.gz",
}


def load_sample_calls(sample_dir, callers):
    """Return {caller: [call dicts]} for one sample directory."""
    calls_by_caller = {}
    for caller in callers:
        pattern = CALLER_FILE_PATTERNS.get(caller, f"*_{caller}.tsv")
        matches = list(Path(sample_dir).rglob(pattern))
        if not matches and pattern.endswith(".vcf.gz"):
            matches = list(Path(sample_dir).rglob(pattern[:-3]))
        if not matches:
            calls_by_caller[caller] = []
            continue
        best = next((m for m in matches if m.parent.name == caller), matches[0])
        calls_by_caller[caller] = CALLER_PARSERS[caller](best)
    return calls_by_caller


# ---------------------------------------------------------------------------
# Candidate regions
# ---------------------------------------------------------------------------

def _merge_candidates(df):
    """Merge overlapping calls (same chr+type) into representative regions.

    After the greedy union merge, calls whose size exceeds MAX_CLUSTER_SIZE_RATIO
    times the median cluster call size are treated as coordinate outliers: they
    are counted in n_callers / callers_supporting but their boundaries are
    excluded from the final candidate coordinates.  This prevents a single
    gene-level caller (e.g. conifer) from inflating a candidate that is
    otherwise supported by exon-level callers at a much smaller scale.
    """
    merged_rows = []
    for (chrom, vtype), group in df.groupby(["chr", "type"]):
        group = group.sort_values("start").reset_index(drop=True)
        used = [False] * len(group)
        for i in range(len(group)):
            if used[i]:
                continue
            caller_i = group.at[i, "caller"]
            cluster = {caller_i}
            cn_by_caller = {caller_i: group.at[i, "cn"] if "cn" in group.columns else None}
            cstart, cend = group.at[i, "start"], group.at[i, "end"]
            intervals = [(group.at[i, "start"], group.at[i, "end"])]
            for j in range(i + 1, len(group)):
                if used[j]:
                    continue
                if overlaps(cstart, cend, group.at[j, "start"],
                            group.at[j, "end"], MIN_OVERLAP_FRAC):
                    caller_j = group.at[j, "caller"]
                    cluster.add(caller_j)
                    if caller_j not in cn_by_caller:
                        cn_by_caller[caller_j] = group.at[j, "cn"] if "cn" in group.columns else None
                    cstart = min(cstart, group.at[j, "start"])
                    cend = max(cend, group.at[j, "end"])
                    intervals.append((group.at[j, "start"], group.at[j, "end"]))
                    used[j] = True
            used[i] = True

            # Clip outlier boundaries when cluster has more than one call.
            if len(intervals) > 1:
                sizes = sorted(e - s for s, e in intervals)
                median_size = sizes[len(sizes) // 2]
                threshold = median_size * MAX_CLUSTER_SIZE_RATIO
                small = [(s, e) for s, e in intervals if e - s <= threshold]
                if small and len(small) < len(intervals):
                    cstart = min(s for s, e in small)
                    cend = max(e for s, e in small)

            cn_label = ";".join(
                "{}:{}".format(c, fmt_cn(cn_by_caller.get(c)))
                for c in sorted(cluster)
            )
            merged_rows.append({"chr": chrom, "start": cstart, "end": cend,
                                "type": vtype, "n_callers": len(cluster),
                                "callers_supporting": ",".join(sorted(cluster)),
                                "cn_label": cn_label})
    return pd.DataFrame(merged_rows)
