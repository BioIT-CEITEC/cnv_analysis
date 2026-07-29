"""
score_samples.py – score caller calls for one or more samples using the trained
BoundedConsensusModel (cnv_consensus_model.pkl).

Expected directory layout:
    sv_dir/
        <sample_name>/
            <caller_files>     (matched by CALLER_FILE_PATTERNS in cnv_consensus_model.py)

Output TSV columns (tab-separated, one row per merged candidate region):
    sample | chr | start | end | type | n_callers | callers_supporting | consensus_score

Usage:
    python score_samples.py \\
        --sv_dir   sv_dir \\
        --model    cnv_consensus_model.pkl \\
        --out      sample_smoothed_variants.tsv \\
        [--cutoff  0.5]

When --cutoff is given a filtered copy is also written alongside --out with
the suffix _filtered_<cutoff>.tsv, which the Nextflow process uses as the
source for _merged_target_consensus.tsv.
"""

import argparse
import pickle
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from cnv_consensus_model import IGNORE_DIRS, load_sample_calls, score_calls  # noqa: E402

EMPTY_HEADER = "sample\tchr\tstart\tend\ttype\tn_callers\tcallers_supporting\tconsensus_score\tcn_label\n"


def load_model(model_path):
    with open(model_path, "rb") as fh:
        return pickle.load(fh)


def main():
    ap = argparse.ArgumentParser(description="Score CNV caller calls with consensus model")
    ap.add_argument("--sv_dir", required=True,
                    help="Directory whose subdirectories are per-sample caller output folders")
    ap.add_argument("--model", required=True,
                    help="Path to the .pkl BoundedConsensusModel file")
    ap.add_argument("--out", required=True,
                    help="Output TSV path for all scored candidates")
    ap.add_argument("--cutoff", type=float, default=None,
                    help="If set, also write a filtered TSV with score >= cutoff")
    args = ap.parse_args()

    model = load_model(args.model)
    callers = model.callers

    sv_dir = Path(args.sv_dir)
    sample_dirs = sorted(
        d for d in sv_dir.iterdir()
        if d.is_dir() and d.name not in IGNORE_DIRS
    )

    frames = []
    for sample_dir in sample_dirs:
        sample = sample_dir.name
        calls_by_caller = load_sample_calls(sample_dir, callers)
        scored = score_calls(calls_by_caller, model, callers)
        if scored.empty:
            continue
        scored.insert(0, "sample", sample)
        frames.append(scored)

    if not frames:
        with open(args.out, "w") as fh:
            fh.write(EMPTY_HEADER)
        return

    out_df = pd.concat(frames, ignore_index=True)
    out_df.to_csv(args.out, sep="\t", index=False)

    if args.cutoff is not None:
        filtered = out_df[out_df["consensus_score"] >= args.cutoff]
        filtered_path = args.out.replace(".tsv", f"_filtered_{args.cutoff:.4f}.tsv")
        filtered.to_csv(filtered_path, sep="\t", index=False)


if __name__ == "__main__":
    main()
