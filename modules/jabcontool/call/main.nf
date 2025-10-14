process JABCONTOOL_CALL {
    publishDir "structural_varcalls/all_samples/jabCoNtool", mode: 'copy'
    conda "${moduleDir}/env.yaml"

    input:
    tuple path(normal_cov), path(normal_snps)
    path organism_regions
    path gc_profile
    path organism_cytoband
    path organism_snps

    output:
    path("calls/final_CNV_probs.tsv"), emit: final_CNV_probs
    path("calls/cohort_info_tab.tsv"), emit: cohort_info_tab

    script:

    def cohort_flag = params.use_cohort_data ? "cohort_data" : "no_previous_cohort_data"
    def snp_bed = params.jabCoNtool_use_snps ? "${organism_snps}" : "no_use_snps"
    def gc_profile_flag = params.jabCoNtool_normalize_to_GC ? "${gc_profile}" : "no_GC_norm"
    def use_cytoband = params.jabCoNtool_remove_centromeres ? "${organism_cytoband}" : "no_cytoband"
    def wgs_or_roi = params.lib_ROI == "wgs" ? "wgs" : "panel"
    //def cov_flag = params.calling_type == "tumor_normal" ? "cov ${tumor_cov} norm_cov ${normal_cov}" : "cov ${normal_cov}"

    """
    mkdir -p calls
    Rscript ${projectDir}/bin/jabConTool_main.R calls/final_CNV_probs.tsv \
        ${organism_regions} \
        ${snp_bed} \
        ${params.calling_type} \
        ${wgs_or_roi} \
        ${gc_profile_flag} \
        ${use_cytoband} \
        ${cohort_flag} \
        ${params.jabCoNtool_predict_TL} \
        ${params.max_CNV_occurance_in_cohort} \
        cov ${normal_cov} 
    """

    stub:
    """
    mkdir -p results
    touch results/final_CNV_probs.tsv
    touch results/cohort_info_tab.tsv
    """
}