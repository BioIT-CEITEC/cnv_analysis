process PREPARE_COHORT_PANELCNMOPS {

    publishDir "cohort_data", mode: 'copy'
    conda "${moduleDir}/env.yaml"

    input:
    tuple path(bam), path(bai)
    path regions_of_interest

    output:
    path("cohort_data/panelcnMOPS_customCohort.RData"), emit: panelcnmops_cohort

    script:
    def bam_files = bam instanceof List ? bam.join(' ') : bam
    """
    mkdir -p cohort_data
    Rscript ${projectDir}/bin/prepare_panelcnMOPS_wrapper.R ${regions_of_interest} cohort_data/panelcnMOPS_customCohort.RData ${bam_files}
    """

    stub:
    """
    mkdir -p cohort_data
    touch cohort_data/panelcnMOPS_customCohort.RData
    """
}