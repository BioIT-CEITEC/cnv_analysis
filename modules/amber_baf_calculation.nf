process AMBER_BAF_CALCULATION {
    tag "$meta.id"
    publishDir "results", mode: 'copy'

    conda (params.conda_enabled ? "bioconda::bioconductor-copynumber" : null)

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("amber/*"), emit: amber

    script:

    def prefix = "${meta.id}"
    
        if (!params.normal_tumor) {
        """
        [ ! -f  ${prefix}.bam ] && ln -s $bams ${prefix}_T.bam
        java -jar $params.tool_dir/amber.jar com.hartwig.hmftools.amber.AmberApplication \
            -tumor ${prefix}_T \
            -tumor_bam ${prefix}_T.bam \
            -output_dir amber/${prefix} \
            -threads $task.cpus \
            -loci $params.organism_amber

        """
    } else {
        """
        [ ! -f  ${prefix}_N.bam ] && ln -s ${bams[0]} ${prefix}_N.bam
        [ ! -f  ${prefix}_T.bam ] && ln -s ${bams[1]} ${prefix}_T.bam
        java -jar $params.tool_dir/amber.jar com.hartwig.hmftools.amber.AmberApplication \
            -reference ${prefix}_N \
            -reference_bam ${prefix}_N.bam \
            -tumor ${prefix}_T \
            -tumor_bam ${prefix}_T.bam \
            -output_dir amber/${prefix} \
            -threads $task.cpus \
            -loci $params.organism_amber
        """
    }

    // TODO: include the threads to the process and the code
}