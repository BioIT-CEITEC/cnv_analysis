include { PREPARE_COHORT_EXOMEDEPTH } from "../../modules/exomeDepth/prepare_cohort/main.nf"
include { CNV_CALL_EXOMEDEPTH } from "../../modules/exomeDepth/cnv_call/main.nf"

workflow EXOMEDEPTH_ANALYSIS {

    take:
    ch_input_bams
    ch_cohort_input
    ch_regions_of_interest
    ch_reference_fasta

    main:

    if (params.use_exomedepth_cohortdata) {
        ch_cohort_data = ch_cohort_input
    } else {
        PREPARE_COHORT_EXOMEDEPTH (
            ch_cohort_input,
            ch_regions_of_interest,
            ch_reference_fasta
        )
        ch_cohort_data = PREPARE_COHORT_EXOMEDEPTH.out.exomeDepth_cohort
    }

    CNV_CALL_EXOMEDEPTH (
        ch_input_bams,
        ch_cohort_data,
        ch_reference_fasta
    )

    emit:
    ch_exomeDepth_varcalls = CNV_CALL_EXOMEDEPTH.out.exomedepth_cnvcalls

}
