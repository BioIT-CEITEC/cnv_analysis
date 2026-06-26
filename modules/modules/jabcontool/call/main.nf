process JABCONTOOL_CALL {
    publishDir "structural_varcalls/all_samples", mode: 'copy'
    conda "${moduleDir}/env.yaml"

    input:
    tuple path(normal_cov), path(normal_snps)
    path organism_regions
    path gc_profile
    path organism_cytoband
    path organism_snps

    output:
    path("jabCoNtool/final_CNV_probs_jabcontool.tsv"), emit: final_CNV_probs
    path("jabCoNtool/cohort_info_tab.tsv"), emit: cohort_info_tab

    script:

    def cohort_flag = params.use_cohort_data ? "cohort_data" : "no_previous_cohort_data"
    def snp_bed = params.jabCoNtool_use_snps ? "${organism_snps}" : "no_use_snps"
    def gc_profile_flag = params.jabCoNtool_normalize_to_GC ? "${gc_profile}" : "no_GC_norm"
    def use_cytoband = params.jabCoNtool_remove_centromeres ? "${organism_cytoband}" : "no_cytoband"
    def wgs_or_roi = params.lib_ROI == "wgs" ? "wgs" : "panel"

    """
    mkdir -p jabCoNtool
    Rscript ${projectDir}/bin/jabConTool_main.R jabCoNtool/final_CNV_probs_jabcontool.tsv \
        ${organism_regions} \
        ${snp_bed} \
        germline \
        ${wgs_or_roi} \
        ${gc_profile_flag} \
        ${use_cytoband} \
        ${cohort_flag} \
        ${params.jabCoNtool_predict_TL} \
        ${params.jabCoNtool_max_CNV_occurance_in_cohort} \
        cov ${normal_cov} 
    """

    stub:
    """
    mkdir -p jabCoNtool
    touch jabCoNtool/final_CNV_probs_jabcontool.tsv
    touch jabCoNtool/cohort_info_tab.tsv
    """
}