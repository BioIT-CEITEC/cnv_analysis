include { PREPARE_COHORT_CNMOPS } from "../../modules/cnMOPS/prepare_cohort/main.nf"
include { CNV_CALL_CNMOPS } from "../../modules/cnMOPS/cnv_call/main.nf"

workflow CNMOPS_ANALYSIS {

    take:
    ch_input_bams
    ch_cohort_bams
    ch_regions_of_interest

    main:

    PREPARE_COHORT_CNMOPS (
        ch_cohort_bams,
        ch_regions_of_interest
    )

ch_cohort_data = PREPARE_COHORT_CNMOPS.out.cnmops_cohort

    CNV_CALL_CNMOPS (
        ch_input_bams,
        ch_cohort_data
    )

}
