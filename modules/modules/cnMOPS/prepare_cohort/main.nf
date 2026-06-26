process PREPARE_COHORT_CNMOPS {

    publishDir "cohort_data", mode: 'copy'
    conda "${moduleDir}/env.yaml"

    input:
    tuple path(bam), path(bai)
    path regions_of_interest

    output:
    path("cnMOPS_customCohort.RData"), emit: cnmops_cohort

    script:

    def bam_files = bam instanceof List ? bam.join(' ') : bam

        """
        Rscript ${projectDir}/bin/prepare_cnMOPS_wrapper.R ${regions_of_interest} cnMOPS_customCohort.RData ${bam_files}
        """

    stub:
        """
        touch cnMOPS_customCohort.RData
        """

}