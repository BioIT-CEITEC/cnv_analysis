process PREPARE_COHORT_PANELCNMOPS {

    publishDir "cohort_data", mode: 'copy'
    conda "${moduleDir}/env.yaml"

    input:
    tuple path(bam), path(bai)
    path regions_of_interest

    output:
    path("panelcnMOPS_customCohort.RData"), emit: panelcnmops_cohort

    script:
    def bam_files = bam instanceof List ? bam.join(' ') : bam
    """
    Rscript ${projectDir}/bin/prepare_panelcnMOPS_wrapper.R ${regions_of_interest} panelcnMOPS_customCohort.RData ${bam_files}
    """

    stub:
    """
    touch panelcnMOPS_customCohort.RData
    """
}