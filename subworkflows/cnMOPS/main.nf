include { PREPARE_COHORT_CNMOPS } from "../../modules/cnMOPS/prepare_cohort/main.nf"
include { CNV_CALL_CNMOPS } from "../../modules/cnMOPS/cnv_call/main.nf"

workflow CNMOPS_ANALYSIS {

    take:
    ch_input_bams
    ch_cohort_input
    ch_regions_of_interest

    main:

    if (params.use_cnmops_cohortdata) {
        ch_cohort_data = ch_cohort_input
    } else {
        PREPARE_COHORT_CNMOPS (
            ch_cohort_input,
            ch_regions_of_interest
        )
        ch_cohort_data = PREPARE_COHORT_CNMOPS.out.cnmops_cohort
    }

    CNV_CALL_CNMOPS (
        ch_input_bams,
        ch_cohort_data
    )

    emit:
    ch_cnmops_varcalls = CNV_CALL_CNMOPS.out.cnmops_cnvcalls

}
