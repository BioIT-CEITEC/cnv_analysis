process MERGE_VARIANT_CALLS {

    publishDir "structural_varcalls/${meta.sample_name}", mode: 'copy'
    conda "${moduleDir}/env.yaml"
    tag "${meta.sample_name}"
    cache false

    input:
    tuple val(meta), path(cnv_varcalls)
    path(dna_panel)
    path(gtf_file)

    output:
    tuple val(meta), path("merged_variants/*_merged_target_consensus.tsv"), path("merged_variants/*_merged_target_consensus.bed"), emit: final_tsvs
    tuple val(meta), path("merged_variants/*_smoothed_variants.tsv"), path("merged_variants/*_smoothed_variants.bed"),       emit: merged_tsv
    tuple val(meta), path("merged_variants/*_per_exon_matrix.tsv"),                                                                emit: raw_matrix

    script:
    def use_jabcontool = params.use_jabcontool ? "--jabcontool_file final_CNV_probs_jabcontool.tsv" : ""
    def use_gatk       = params.use_gatk       ? "--gatk_vcf_file ${meta.sample_name}_gatk.vcf.gz"  : ""
    """
    Rscript ${projectDir}/bin/merging_and_smoothing.R \
        --input_dir    "./" \
        --capture_bed  ${dna_panel} \
        --gtf_file     ${gtf_file} \
        --output_dir   "merged_variants" \
        --sample_id    ${meta.sample_name} \
        --min_callers  ${params.min_number_callers} \
        --max_call_size ${params.size_threshold} \
        --smooth_gap_bp ${params.smooth_gap_bp} \
        ${use_jabcontool} \
        ${use_gatk}
    """

    stub:
        """
        mkdir -p merged_variants
        touch merged_variants/${meta.sample_name}_merged_target_consensus.tsv
        touch merged_variants/${meta.sample_name}_merged_target_consensus.bed
        touch merged_variants/${meta.sample_name}_smoothed_variants.tsv
        touch merged_variants/${meta.sample_name}_smoothed_variants.bed
        touch merged_variants/${meta.sample_name}_per_exon_matrix.tsv
        """

}