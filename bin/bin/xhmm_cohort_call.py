import argparse
import os
import re
import shlex
import subprocess
import tempfile


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--bams",         nargs="+", required=True)
    p.add_argument("--sample-names", required=True,
                   help="Comma-separated list of sample names, in the same order as --bams")
    p.add_argument("--bed",          required=True)
    p.add_argument("--output",       required=True)
    p.add_argument("--threads",      type=int, default=1)
    p.add_argument("--log",          default="xhmm_call.log")
    # Filtering params
    p.add_argument("--min-target-size",   type=int,   default=10)
    p.add_argument("--max-target-size",   type=int,   default=10000)
    p.add_argument("--min-mean-target-rd",type=float, default=10.0)
    p.add_argument("--max-mean-target-rd",type=float, default=500.0)
    p.add_argument("--min-mean-sample-rd",type=float, default=25.0)
    p.add_argument("--max-mean-sample-rd",type=float, default=200.0)
    p.add_argument("--max-sd-sample-rd",  type=float, default=150.0)
    p.add_argument("--pve-mean-factor",   type=float, default=0.7)
    p.add_argument("--max-sd-target-rd",  type=float, default=30.0)
    p.add_argument("--discover-some-qual-thresh", type=float, default=29.99)
    p.add_argument("--params-values",  nargs=9, type=float,
                   default=[1e-8, 6, 70, -3.0, 1.0, 0.0, 1.0, 3.0, 1.0])
    return p.parse_args()


def _q(x):
    return shlex.quote(str(x))


def _run(cmd, log_filename):
    with open(log_filename, "at") as lf:
        lf.write(f"\n## CMD: {cmd}\n")
    subprocess.run(cmd + " >> " + _q(log_filename) + " 2>&1", shell=True, check=True)


def _matrix_dimensions(matrix_file):
    with open(matrix_file) as f:
        header = f.readline().strip().split("\t")
        if len(header) < 2:
            return 0, 0
        num_targets = len(header) - 1
        num_samples = sum(1 for line in f if line.strip())
    return num_samples, num_targets


def _parse_mean_coverage(line):
    fields = line.strip().split()
    if len(fields) < 6:
        return None
    mean_token = fields[-2] if len(fields) >= 7 else fields[-1]
    m = re.match(r"^([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)", mean_token)
    return m.group(1) if m else None


def main():
    args = parse_args()

    log_filename  = os.path.abspath(args.log)
    bed_file      = os.path.abspath(args.bed)
    bam_files     = [os.path.abspath(p) for p in args.bams]
    sample_names  = args.sample_names.split(",")

    with open(log_filename, "wt") as f:
        f.write("\n##\n## XHMM cohort CNV call\n##\n")

    if len(bam_files) != len(sample_names):
        raise ValueError(f"--bams ({len(bam_files)}) and --sample-names ({len(sample_names)}) have different lengths")

    xcnv_out = os.path.abspath(args.output)
    out_dir  = os.path.dirname(xcnv_out)
    work_dir = os.path.join(out_dir, "xhmm_work")
    os.makedirs(work_dir, exist_ok=True)
    os.makedirs(out_dir,  exist_ok=True)

    # Build XHMM-format intervals from BED
    intervals = []
    with open(bed_file) as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                continue
            chrom = fields[0]
            start = int(fields[1]) + 1
            end   = int(fields[2])
            if end < start:
                continue
            intervals.append(f"{chrom}:{start}-{end}")

    if not intervals:
        raise ValueError("No valid intervals found in BED file")

    # Per-sample sambamba coverage → read-depth matrix
    rd_matrix    = os.path.join(work_dir, "DATA.RD.txt")
    sample_depths = {}

    for sample, bam in zip(sample_names, bam_files):
        if not os.path.exists(bam + ".bai"):
            _run("samtools index " + _q(bam), log_filename)

        tmp_cov = tempfile.NamedTemporaryFile(
            prefix="xhmm_cov_", suffix=".txt", dir=work_dir, delete=False
        )
        tmp_cov.close()

        subprocess.run(
            "sambamba depth region --regions=" + _q(bed_file) +
            " --min-base-quality=0 --nthreads=" + _q(args.threads) +
            " " + _q(bam) + " > " + _q(tmp_cov.name) + " 2>&1",
            shell=True, check=True
        )

        values      = []
        total_lines = 0
        first_lines = []

        with open(tmp_cov.name) as cov_f:
            for raw in cov_f:
                total_lines += 1
                if len(first_lines) < 5:
                    first_lines.append(raw.rstrip("\n"))
                if not raw or raw.startswith("#"):
                    continue
                mean_cov = _parse_mean_coverage(raw)
                if mean_cov is not None:
                    values.append(mean_cov)

        if len(values) != len(intervals):
            with open(log_filename, "at") as f:
                f.write(
                    f"XHMM debug: sample={sample} total_lines={total_lines} "
                    f"parsed={len(values)} intervals={len(intervals)} cov_file={tmp_cov.name}\n"
                )
                for line in first_lines:
                    f.write(f"XHMM debug line: {line}\n")
            raise ValueError(
                f"Coverage row count mismatch for '{sample}' ({len(values)} vs {len(intervals)} intervals)"
            )

        os.remove(tmp_cov.name)
        sample_depths[sample] = values

    with open(rd_matrix, "w") as f:
        f.write("SAMPLE\t" + "\t".join(intervals) + "\n")
        for sample in sample_names:
            f.write(sample + "\t" + "\t".join(sample_depths[sample]) + "\n")

    # XHMM pipeline
    filtered_centered  = os.path.join(work_dir, "DATA.filtered_centered.RD.txt")
    filtered_targets_1 = filtered_centered + ".filtered_targets.txt"
    filtered_samples_1 = filtered_centered + ".filtered_samples.txt"

    pca_base       = os.path.join(work_dir, "DATA.RD_PCA")
    pca_normalized = os.path.join(work_dir, "DATA.PCA_normalized.txt")

    zscores            = os.path.join(work_dir, "DATA.PCA_normalized.filtered.sample_zscores.RD.txt")
    filtered_targets_2 = zscores + ".filtered_targets.txt"
    filtered_samples_2 = zscores + ".filtered_samples.txt"

    same_filtered    = os.path.join(work_dir, "DATA.same_filtered.RD.txt")
    aux_xcnv         = os.path.join(work_dir, "DATA.aux_xcnv")
    posterior_prefix = os.path.join(work_dir, "DATA.posterior")

    params_file = os.path.join(work_dir, "params.txt")
    with open(params_file, "w") as pf:
        pf.write("\t".join(str(v) for v in args.params_values) + "\n")

    def _filter_matrix(out_file, out_targets, out_samples, extra_flags=""):
        _run(
            "xhmm --matrix"
            " -r " + _q(rd_matrix) +
            " --centerData --centerType target"
            " -o " + _q(out_file) +
            " --outputExcludedTargets " + _q(out_targets) +
            " --outputExcludedSamples " + _q(out_samples) +
            extra_flags,
            log_filename
        )

    _filter_matrix(
        filtered_centered, filtered_targets_1, filtered_samples_1,
        " --minTargetSize "    + _q(args.min_target_size) +
        " --maxTargetSize "    + _q(args.max_target_size) +
        " --minMeanTargetRD "  + _q(args.min_mean_target_rd) +
        " --maxMeanTargetRD "  + _q(args.max_mean_target_rd) +
        " --minMeanSampleRD "  + _q(args.min_mean_sample_rd) +
        " --maxMeanSampleRD "  + _q(args.max_mean_sample_rd) +
        " --maxSdSampleRD "    + _q(args.max_sd_sample_rd)
    )

    n_samples, n_targets = _matrix_dimensions(filtered_centered)
    if n_targets <= n_samples:
        with open(log_filename, "at") as f:
            f.write(
                f"XHMM warning: strict filtering gave {n_targets} targets for {n_samples} samples; "
                "retrying with relaxed filters.\n"
            )
        _filter_matrix(
            filtered_centered, filtered_targets_1, filtered_samples_1,
            " --minTargetSize 1 --maxTargetSize 1000000000"
            " --minMeanTargetRD 0 --maxMeanTargetRD 1000000000"
            " --minMeanSampleRD 0 --maxMeanSampleRD 1000000000"
            " --maxSdSampleRD 1000000000"
        )
        n_samples, n_targets = _matrix_dimensions(filtered_centered)
        if n_targets <= n_samples:
            raise ValueError(
                f"Insufficient targets after adaptive filtering ({n_targets} targets, {n_samples} samples)"
            )

    _run("xhmm --PCA -r " + _q(filtered_centered) + " --PCAfiles " + _q(pca_base), log_filename)

    _run(
        "xhmm --normalize"
        " -r " + _q(filtered_centered) +
        " --PCAfiles " + _q(pca_base) +
        " --normalizeOutput " + _q(pca_normalized) +
        " --PCnormalizeMethod PVE_mean"
        " --PVE_mean_factor " + _q(args.pve_mean_factor),
        log_filename
    )

    _run(
        "xhmm --matrix"
        " -r " + _q(pca_normalized) +
        " --centerData --centerType sample --zScoreData"
        " -o " + _q(zscores) +
        " --outputExcludedTargets " + _q(filtered_targets_2) +
        " --outputExcludedSamples " + _q(filtered_samples_2) +
        " --maxSdTargetRD " + _q(args.max_sd_target_rd),
        log_filename
    )

    _run(
        "xhmm --matrix"
        " -r " + _q(rd_matrix) +
        " --excludeTargets " + _q(filtered_targets_1) +
        " --excludeTargets " + _q(filtered_targets_2) +
        " --excludeSamples " + _q(filtered_samples_1) +
        " --excludeSamples " + _q(filtered_samples_2) +
        " -o " + _q(same_filtered),
        log_filename
    )

    _run(
        "xhmm --discover"
        " -p " + _q(params_file) +
        " -r " + _q(zscores) +
        " -R " + _q(same_filtered) +
        " -c " + _q(xcnv_out) +
        " -a " + _q(aux_xcnv) +
        " -s " + _q(posterior_prefix) +
        " -t " + _q(args.discover_some_qual_thresh),
        log_filename
    )


if __name__ == "__main__":
    main()
