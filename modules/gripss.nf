process GRIPSS_CALL {
    tag "$meta.id"
    publishDir "results/gridss/$meta.id", mode: 'copy'

    conda (params.conda_enabled ? "bioconda::bioconductor-copynumber" : null)

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.somatic.vcf.gz"), path("*.sv.vcf.gz") , emit: vcfs

    script:

    def prefix = "${meta.id}"

        if (!params.normal_tumor) {
        """
        java -jar $params.tool_dir/gripss.jar \
            -sample ${prefix}_T \
            -ref_genome_version $params.organism_genome_version \
            -ref_genome $params.organism_fasta \
            -pon_sgl_file gridss_pon_single_breakend.bed \
            -pon_sv_file gridss_pon_breakpoint.bedpe \
            -vcf ${prefix}.gridss.unfiltered.vcf.gz

        mv ${prefix}_T.somatic.filtered.vcf.gz ${prefix}_T.sv.vcf.gz

        """
    } else {
        """
        java -jar $params.tool_dir/gripss.jar \
            -sample ${prefix}_T \
            -reference ${prefix}_N \
            -ref_genome_version $params.organism_genome_version \
            -ref_genome $params.organism_fasta \
            -pon_sgl_file gridss_pon_single_breakend.bed \
            -pon_sv_file gridss_pon_breakpoint.bedpe \
            -vcf ${prefix}.gridss.unfiltered.vcf.gz

        mv ${prefix}_T.gripss.somatic.filtered.vcf.gz ${prefix}_T.sv.vcf.gz

        """
    }

}