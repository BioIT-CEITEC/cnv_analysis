process MERGE_VARIANT_CALLS {

    publishDir "structural_varcalls/${meta.sample_name}", mode: 'copy'
    conda "${moduleDir}/env.yaml"
    tag "${meta.sample_name}"

    input:
    tuple val(meta), path(cnv_varcalls)
    path(dna_panel)

    output:
    tuple val(meta), path("merged_variants/*_merged_target_consensus.tsv"), path("merged_variants/*.bed"), emit: merged_calls
    tuple val(meta), path("merged_variants/*_smoothed_variants.tsv"),                              emit: smoothed_calls
    tuple val(meta), path("merged_variants/*_raw_callers_matrix.tsv"),                             emit: raw_matrix

    script:

        """
        Rscript ${projectDir}/bin/merging_wrapper.R ${dna_panel} ${meta.sample_name} ${params.min_number_callers}
        """

    stub:
        """
        mkdir -p merged_variants
        touch merged_variants/target.bed
        touch merged_variants/antitarget.bed
        """

}