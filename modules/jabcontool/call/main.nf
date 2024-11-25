process JABCONTOOL_CALL {

    input:
    tuple path(tumor_cov), path(tumor_snps), path(normal_cov), path(normal_snps)
    path organism_dna_panel
    path organism_snps
    path gc_profile
    path binned_genome
    path organism_cytoband
    path cohort_data

    output:
    path("final_CNV_probs.tsv"), emit: final_CNV_probs
    path("cohort_info_tab.tsv"), emit: cohort_info_tab

    script:

    def cohort_flag = params.use_cohort_data ? "cohort_info_tab.tsv" : "no_previous_cohort_data"
    def region_bed = params.lib_ROI == "wgs" ? "${binned_genome}" : "${organism_dna_panel}"
    def snp_bed = params.jabCoNtool_use_snps ? "${organism_snps}" : "no_use_snps"
    def gc_profile_flag = params.jabCoNtool_normalize_to_GC ? "${gc_profile}" : "no_GC_norm"
    def use_cytoband = params.jabCoNtool_remove_centromeres ? "${organism_cytoband}" : "no_cytoband"
    def wgs_or_roi = params.lib_ROI == "wgs" ? "wgs" : "panel"
    def normal_cov_flag = params.calling_type == "tumor_normal" ? "norm_cov ${normal_cov}" : ""

    """
    Rscript jabConTool_main.R results/final_CNV_probs.tsv \
        ${region_bed} \
        ${snp_bed} \
        ${params.calling_type} \
        ${wgs_or_roi} \
        ${gc_profile_flag} \
        ${use_cytoband} \
        ${cohort_flag} \
        ${params.jabCoNtool_predict_TL} \
        ${params.max_CNV_occurance_in_cohort} \
        cov ${tumor_cov} \
        ${normal_cov_flag}
    """

    stub:
    """
    mkdir -p results
    touch results/final_CNV_probs.tsv
    touch results/cohort_info_tab.tsv
    """
}