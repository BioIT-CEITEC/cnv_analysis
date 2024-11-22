process COVERAGE_CALC {

    input:
    tuple val(meta), path(bams)
    path(binned_genome), path(gc_profile)

    output:
    tuple val(meta), path("*.region_coverage.tsv"), emit: region_coverage

    script:

    def prefix = "${meta.id}"
    def region_bed = params.normal_tumor ? "binned_genome_${params.window_size}.bed" : "${params.organism_dna_panel}"

        if (!params.normal_tumor) {
        """
        [ ! -f  ${prefix}.bam ] && ln -s $bams ${prefix}_T.bam
        bedtools coverage -sorted -a $region_bed -b ${prefix}_T.bam -g $params.organism_dict > ${prefix}_T.region_coverage.tsv

        """
    } else {
        """
        [ ! -f  ${prefix}_N.bam ] && ln -s ${bams[0]} ${prefix}_N.bam
        [ ! -f  ${prefix}_T.bam ] && ln -s ${bams[1]} ${prefix}_T.bam
        bedtools coverage -sorted -a $region_bed -b ${prefix}_T.bam -g $params.organism_dict > ${prefix}_T.region_coverage.tsv
        bedtools coverage -sorted -a $region_bed -b ${prefix}_N.bam -g $params.organism_dict > ${prefix}_N.region_coverage.tsv
        """
    }

    // TODO: include the threads to the process and the code
}

process SNP_AF_CALC {

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("*_T.snpAF.tsv"), path("*_N.snpAF.tsv"), emit: snpAF

    script:

    def prefix = "${meta.id}"

        if (!params.normal_tumor) {
        """
        [ ! -f  ${prefix}.bam ] && ln -s $bams ${prefix}_T.bam

        alleleCounter -r $params.organism_fasta -l $params.organism_snps_tsv -b ${prefix}_T.bam -o ${prefix}_T.snpAF.tsv
        touch dummy_N.snpAF.tsv
        """
    } else {
        """
        [ ! -f  ${prefix}_N.bam ] && ln -s ${bams[0]} ${prefix}_N.bam
        [ ! -f  ${prefix}_T.bam ] && ln -s ${bams[1]} ${prefix}_T.bam

        alleleCounter -r $params.organism_fasta -l $params.organism_snps_tsv -b ${prefix}_N.bam -o ${prefix}_N.snpAF.tsv
        alleleCounter -r $params.organism_fasta -l $params.organism_snps_tsv -b ${prefix}_T.bam -o ${prefix}_T.snpAF.tsv
        """
    }

    // TODO: include the threads to the process and the code
}

process JABCONTOOL_CALL {

    input:
    path(normal)
    path(tumor)
    path(cohort_data)

    output:

    script:

    def cohort_data = params.use_cohort_data ? "cohort_info_tab.tsv" : "no_previous_cohort_data"
    def region_bed = params.lib_ROI == "wgs" ? "binned_genome_${params.window_size}.bed" : "${params.organism_dna_panel}"
    def snp_bed = params.jabCoNtool_use_snps ? "${params.organism_snps_tsv}" : "no_use_snps"
    def gc_profile = params.jabCoNtool_normalize_to_GC ? "GC_profile_${params.windows_size}.cnp" : "no_GC_norm"
    def cytoband = params.jabCoNtool_remove_centromeres ? "${params.organism_cytoband}" : "no_cytoband"
    def wgs_or_roi = params.lib_ROI == "wgs" ? "wgs" : "panel"
    def normal_cov = params.calling_type == "tumor_normal" ? "norm_cov ${normal}" : ""

    """
    Rscript jabConTool_main.R final_CNV_probs.tsv \
        $region_bed \
        $snp_bed \
        $params.calling_type \
        $wgs_or_roi \
        $gc_profile \
        $cytoband \
        $cohort_data \
        $params.jabCoNtool_predict_TL \
        $params.max_CNV_occurance_in_cohort \
        cov $tumor \
        $normal_cov
    """
}