process CNV_CALL_FREEC {
    tag "${meta.sample_name}"
    conda "${moduleDir}/../env.yaml"

    input:
    tuple val(meta), path(bam), path(bam_bai)
    path chrlen
    path mappability_bg
    path capture_bed

    output:
    tuple val(meta),
          path("${meta.sample_name}_freec/*.bam_CNVs"),
          path("${meta.sample_name}_freec/${meta.sample_name}_freec_config.txt"),
          path("${meta.sample_name}_freec/*.bam_ratio.txt"),
          emit: var_call

    script:
    def window_line = (params.freec_window && params.freec_window != "0") \
        ? "config_lines.append('window = ${params.freec_window}')" : ""

    """
    mkdir -p ${meta.sample_name}_freec

    python3 << CODE
import os

bed_chroms = set()
with open("${capture_bed}") as fh:
    for line in fh:
        if not line or line.startswith("#"):
            continue
        fields = line.rstrip("\\n").split("\\t")
        if fields and fields[0]:
            bed_chroms.add(fields[0])

filtered_chrlen = "${meta.sample_name}_freec/${meta.sample_name}_chrLen.filtered.txt"
kept = 0
with open("${chrlen}") as src, open(filtered_chrlen, "w") as dst:
    for line in src:
        fields = line.rstrip("\\n").split("\\t")
        if len(fields) >= 2 and fields[0] in bed_chroms:
            dst.write(line)
            kept += 1

if kept == 0:
    raise ValueError(
        "No overlapping chromosomes between chrLenFile and capture BED for sample ${meta.sample_name}"
    )

config_lines = [
    "[general]",
    "chrLenFile = "               + filtered_chrlen,
    "ploidy = ${params.freec_ploidy}",
    "outputDir = ${meta.sample_name}_freec/",
    "gemMappabilityFile = ${mappability_bg}",
    "minMappabilityPerWindow = ${params.freec_min_map}",
    "minExpectedGC = ${params.freec_min_expected_gc}",
    "maxExpectedGC = ${params.freec_max_expected_gc}",
    "coefficientOfVariation = ${params.freec_coeff_var}",
    "breakPointThreshold = ${params.freec_breakpoint_threshold}",
    "breakPointType = 4",
    "maxThreads = ${task.cpus}",
    "forceGCcontentNormalization = 1",
]
${window_line}

config_lines += [
    "",
    "[sample]",
    "mateFile = ${bam}",
    "inputFormat = BAM",
    "mateOrientation = 0",
    "",
    "[target]",
    "captureRegions = ${capture_bed}",
]

config_file = "${meta.sample_name}_freec/${meta.sample_name}_freec_config.txt"
with open(config_file, "w") as fh:
    fh.write("\\n".join(config_lines) + "\\n")
CODE

    freec -conf ${meta.sample_name}_freec/${meta.sample_name}_freec_config.txt
    """

    stub:
    """
    touch ${meta.sample_name}_freec/${meta.sample_name}_${meta.sample_name}.bam_CNVs
    touch ${meta.sample_name}_freec/${meta.sample_name}_freec_config.txt
    touch ${meta.sample_name}_freec/${meta.sample_name}.bam_ratio.txt
    """

}
