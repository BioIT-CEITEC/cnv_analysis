include { PREPARE_COHORT_PANELCNMOPS } from "../../modules/panelcn.MOPS/prepare_cohort/main.nf"
include { CNV_CALL_PANELCNMOPS } from "../../modules/panelcn.MOPS/cnv_call/main.nf"

workflow PANELCNMOPS_ANALYSIS {

    take:
    ch_input_bams
    ch_cohort_input
    ch_regions_of_interest

    main:

    if (params.use_panelcnmops_cohortdata) {
        ch_cohort_data = ch_cohort_input
    } else {
        PREPARE_COHORT_PANELCNMOPS (
            ch_cohort_input,
            ch_regions_of_interest
        )
        ch_cohort_data = PREPARE_COHORT_PANELCNMOPS.out.panelcnmops_cohort
    }

    CNV_CALL_PANELCNMOPS (
        ch_input_bams,
        ch_cohort_data
    )

    emit:
    ch_panelcnmops_varcalls = CNV_CALL_PANELCNMOPS.out.panelcnMOPS_cnvcalls

}
