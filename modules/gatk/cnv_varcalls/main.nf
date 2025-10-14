process CNV_VARCALLS_GATK {

    publishDir "structural_varcalls/all_samples/gatk", mode: 'copy'
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
        mkdir -p cohort_data

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
        mkdir -p cohort_data
        touch cohort_data/target.bed
        touch cohort_data/antitarget.bed
        """

}