process AMBER_BAF_CALCULATION {
    tag "${meta.id}"
    publishDir "results/purple/amber/${meta.id}", mode: 'copy'

    conda (params.conda_enabled ? "bioconda::bioconductor-copynumber" : null)

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    val genome_version
    path heterozygous_sites
    path lib_ROI
    output:
    tuple val(meta), path("amber/"), emit: amber

    script:

    def normal_arg = params.normal_tumor ? "-reference \"${meta.id}_N\"" : ""

    def normal_bam = params.normal_tumor ? "-reference_bam \"${meta.id}_N.bam\"" : ""

    def regions_of_interest = (params.lib_ROI && params.lib_ROI != "wgs") ? "-target_regions_bed ${ROI_bed}" : ""

        if (!params.normal_tumor) {
        """
        [ ! -f  ${prefix}.bam ] && ln -s $bams ${prefix}_T.bam
        java -jar $params.tool_dir/cnv_tools/amber.jar com.hartwig.hmftools.amber.AmberApplication \
            -tumor ${prefix}_T \
            -tumor_bam ${prefix}_T.bam \
            -output_dir amber/ \
            -threads $task.cpus \
            -loci $params.organism_amber

        """
    } else {
        """
        [ ! -f  ${prefix}_N.bam ] && ln -s ${bams[0]} ${prefix}_N.bam
        [ ! -f  ${prefix}_T.bam ] && ln -s ${bams[1]} ${prefix}_T.bam
        java -jar $params.tool_dir/cnv_tools/amber.jar com.hartwig.hmftools.amber.AmberApplication \
            -reference ${prefix}_N \
            -reference_bam ${prefix}_N.bam \
            -tumor ${prefix}_T \
            -tumor_bam ${prefix}_T.bam \
            -output_dir amber/ \
            -threads $task.cpus \
            -loci $params.organism_amber
        """
    }

    // TODO: include the threads to the process and the code
}