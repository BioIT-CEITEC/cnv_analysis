process PURPLE_CALL {
    tag "$meta.id"
    publishDir "results/purple/${meta.id}", mode: 'copy'

    conda (params.conda_enabled ? "bioconda::bioconductor-copynumber" : null)

    input:
    tuple val(meta), path(amber), path(cobalt), path(somatic), path(sv)

    output:
    tuple val(meta), path("*.purity.tsv"), path("*.purity.qc"), emit: purity
    tuple val(meta), path("*.cnv.somatic.tsv"), path("*.cnv.gene.tsv"), emit: cnvs
    tuple val(meta), path("*.sv.vcf.gz"), path("*.somatic.vcf.gz"), emit: vcfs



    script:

    def prefix = "${meta.id}"
    def amber_folder = amber[0].getParent()
    def cobalt_folder = cobalt[0].getParent()
    
        if (!params.normal_tumor) {
        """
        java -jar $params.tool_dir/cnv_tools/purple.jar \
            -tumor ${prefix}_T  \
            -amber $amber_folder \
            -cobalt $cobalt_folder \
            -gc_profile $params.organism_cobalt \
            -ref_genome_version $params.organism_genome_version \
            -ref_genome $params.organism_fasta \
            -ensembl_data_dir $params.organism_ensembl_cache \
            -somatic_vcf ${prefix}_T.gripss.somatic.vcf.gz \
            -structural_vcf ${prefix}_T.gripss.sv.vcf.gz \
            -circos $params.tool_dir/circos-0.69-6/bin/circos 

        """
    } else {
        """
        java -jar $params.tool_dir/cnv_tools/purple.jar \
            -reference ${prefix}_N  \
            -tumor ${prefix}_T  \
            -amber $amber_folder \
            -cobalt $cobalt_folder \
            -gc_profile $params.organism_cobalt \
            -ref_genome_version $params.organism_genome_version \
            -ref_genome $params.organism_fasta \
            -ensembl_data_dir $params.organism_ensembl_cache \
            -somatic_vcf ${prefix}_T.gripss.somatic.vcf.gz \
            -structural_vcf ${prefix}_T.gripss.sv.vcf.gz \
            -circos $params.tool_dir/circos-0.69-6/bin/circos 
        """
    }

}