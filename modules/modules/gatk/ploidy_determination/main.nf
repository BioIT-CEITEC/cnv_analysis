process GERMLINE_PLOIDY_DETERMINATION_GATK {

    conda "${moduleDir}/../env.yaml"

    input:
    path qc_filtered_intervals
    path read_counts
    path ploidy_priors
 
    output:
    path("model/*"), emit: ploidy_model_gatk

    script:
    """
    mkdir -p model
    gatk --java-options "-Xmx10g" DetermineGermlineContigPloidy \\
        -L ${qc_filtered_intervals} \\
        --interval-merging-rule OVERLAPPING_ONLY \\
        -I ${read_counts.join(" -I ")} \\
        --contig-ploidy-priors ${ploidy_priors} \\
        --output model/ \\
        --output-prefix ploidy \\
        --verbosity DEBUG
    """
    stub:
    """
    mkdir -p model
    touch model/placeholder.txt
    """
}













