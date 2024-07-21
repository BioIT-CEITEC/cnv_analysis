process COBALT_READ_DEPTH {
    tag "$meta.id"
    publishDir "results/purple/cobalt", mode: 'copy'

    conda (params.conda_enabled ? "bioconda::bioconductor-copynumber" : null)

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("cobalt/*"), emit: tsv

    script:

    def prefix = "${meta.id}"
    
        if (!params.normal_tumor) {
        """
        [ ! -f  ${prefix}.bam ] && ln -s $bams ${prefix}_T.bam
        java -jar -Xmx8G $params.tool_dir/cobalt.jar \
            -tumor ${prefix}_T \
            -tumor_bam ${prefix}_T.bam \
            -output_dir cobalt/${prefix} \
            -threads $task.cpus \
            -gc_profile $params.organism_cobalt
            -tumor-only-diploid-bed $params.organism_cobalt_tumor

        """
    } else {
        """
        [ ! -f  ${prefix}_N.bam ] && ln -s ${bams[0]} ${prefix}_N.bam
        [ ! -f  ${prefix}_T.bam ] && ln -s ${bams[1]} ${prefix}_T.bam
        java -jar -Xmx8G $params.tool_dir/cobalt.jar \
            -reference ${prefix}_N \
            -reference_bam ${prefix}_N.bam \
            -tumor ${prefix}_T \
            -tumor_bam ${prefix}_T.bam \
            -output_dir amber/${prefix} \
            -threads $task.cpus \
            -gc_profile $params.organism_cobalt
        """
    }

    // TODO: include the threads to the process and the code
}