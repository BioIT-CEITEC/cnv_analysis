process COHORT_GENERATION_GATK {

    conda "${moduleDir}/../env.yaml"

    input:
    path interval_list
    path read_counts
    path annotated_intervals
    path ploidy_calls

    output:
    path("cohort_data/*"), emit: cohort_model

    script:

    """
    mkdir -p cohort_data

    gatk --java-options "-Xmx20g" GermlineCNVCaller \\
        --run-mode COHORT \\
        -L ${interval_list} \\
        -I ${read_counts.join(" -I ")} \\
        --contig-ploidy-calls ploidy-calls/ \\
        --annotated-intervals ${annotated_intervals} \\
        --interval-merging-rule OVERLAPPING_ONLY \\
        --output cohort_data \\
        --output-prefix cohort \\
        --verbosity DEBUG
    """

    stub:
    """
    mkdir -p cohort_data
    touch cohort_data/target.bed
    """

}