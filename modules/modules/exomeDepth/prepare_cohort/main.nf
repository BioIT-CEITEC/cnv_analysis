process PREPARE_COHORT_EXOMEDEPTH {

    publishDir "cohort_data", mode: 'copy'
    conda "${moduleDir}/env.yaml"

    input:
    tuple path(bam), path(bai)
    path regions_of_interest
    path reference_fasta

    output:
    path("exomeDepth_customCohort.RData"), emit: exomeDepth_cohort

    script:

    def bam_files = bam instanceof List ? bam.join(' ') : bam

    """
    Rscript ${projectDir}/bin/prepare_exomeDepth_wrapper.R ${regions_of_interest} ${reference_fasta} exomeDepth_customCohort.RData ${bam_files} 
    """

    stub:
        """
        touch exomeDepth_customCohort.RData
        """

}