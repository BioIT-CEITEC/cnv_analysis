process GRIDSS_CALL {
    tag "$meta.id"
    publishDir [ { file -> file.name.endsWith('.unfiltered.vcf.gz') ? "$meta.id" : "PON" } ], mode: 'copy'

    conda (params.conda_enabled ? "bioconda::bioconductor-copynumber" : null)

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("*.unfiltered.vcf.gz"), emit: vcf
    tuple val(meta), path("gridss_pon_breakpoint.bedpe"), path("gridss_pon_single_breakend.bed"), emit: bedpe

    script:

    def prefix = "${meta.id}"

        if (!params.normal_tumor) {
        """
        [ ! -f  ${prefix}.bam ] && ln -s $bams ${prefix}_T.bam
        gridss --jvmheap 30g \
            --reference $params.organism_fasta \
            --output ${prefix}.gridss.unfiltered.vcf.gz \
            --threads $task.cpus \
            --labels ${prefix}_T \
            ${prefix}_T.bam

        java -Xmx8G \
            -cp $params.tool_dir/cnv_tools/gridss-2.13.1-gridss-jar-with-dependencies.jar
            gridss.GeneratePonBedpe \
            $(ls -1 *unfiltered.vcf.gz | awk ' { print "INPUT=" $0 }' | head -$n) \
            O=gridss_pon_breakpoint.bedpe \
            SBO=gridss_pon_single_breakend.bed \
            REFERENCE=$params.organism_fasta

        """
    } else {
        """
        [ ! -f  ${prefix}_N.bam ] && ln -s ${bams[0]} ${prefix}_N.bam
        [ ! -f  ${prefix}_T.bam ] && ln -s ${bams[1]} ${prefix}_T.bam
        gridss --jvmheap 30g \
            --reference $params.organism_fasta \
            --output ${prefix}.gridss.unfiltered.vcf.gz \
            --threads $task.cpus \
            --labels ${prefix}_N,${prefix}_T \
            ${prefix}_N.bam ${prefix}_T.bam

        java -Xmx8G \
            -cp $params.tool_dir/cnv_tools/gridss-2.13.1-gridss-jar-with-dependencies.jar
            gridss.GeneratePonBedpe \
            $(ls -1 *unfiltered.vcf.gz | awk ' { print "INPUT=" $0 }' | head -$n) \
            O=gridss_pon_breakpoint.bedpe \
            SBO=gridss_pon_single_breakend.bed \
            REFERENCE=$params.organism_fasta

        """
    }

}