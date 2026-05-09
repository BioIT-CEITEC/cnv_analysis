process CNV_VARCALLS_GATK {
    publishDir "structural_varcalls/${meta.sample_name}/gatk", mode: 'copy'
    conda "${moduleDir}/../env.yaml"
    tag "${meta.sample_name}"

    input:
    tuple val(meta), val(read_counts)
    path ploidy_calls
    path cohort_data

    output:
    tuple val(meta), path("germline_calls/*"), emit: germline_calls

    script:
    """
    gatk GermlineCNVCaller \\
        -I ${read_counts} \\
        --run-mode CASE \\
        --contig-ploidy-calls ploidy-calls/ \\
        --model cohort-model \\
        --output germline_calls \\
        --output-prefix germline \\
        --verbosity DEBUG
    """

    stub:
    """
    mkdir -p germline_calls
    touch germline_calls/placeholder.vcf.gz
    """

}