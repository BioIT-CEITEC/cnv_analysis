process MERGE_VARIANT_CALLS {

    publishDir "merged_varcalls/", mode: 'copy'
    conda "${moduleDir}/env.yaml"
    tag "${meta.sample_name}"

    input:
    tuple val(meta), path(cnv_varcalls)
    path(dna_panel)

    output:
    tuple val(meta), path("results/*.tsv"), path("results/*.bed"), emit: merged_calls

    script:

        """
        Rscript ${projectDir}/bin/merging_wrapper.R ${dna_panel} ${meta.sample_name} ${params.min_number_callers}
        """

    stub:
        """
        mkdir -p cohort_data
        touch cohort_data/target.bed
        touch cohort_data/antitarget.bed
        """

}