include { CNV_CALL_XHMM } from "../../modules/XHMM/cnv_call/main.nf"
include { PER_SAMPLE_CALL_XHMM } from "../../modules/XHMM/per_sample_call/main.nf"

workflow XHMM_ANALYSIS {

    take:
    ch_input_bams
    ch_regions_of_interest

    main:

    ch_input_bams
        .map { sample ->
            sample[1..-1]
        }
        .flatten()
        .collect()
        .set { ch_all_bam_files }

    ch_input_bams
        .map { meta, bam, bai -> bam.simpleName }
        .collect()
        .set { ch_unique_basenames }

    ch_unique_basenames.view()

    CNV_CALL_XHMM (
        ch_all_bam_files,
        ch_regions_of_interest,
        ch_unique_basenames
    )

    ch_cohort_varcalls = CNV_CALL_XHMM.out.xhmm_cnvcalls

    PER_SAMPLE_CALL_XHMM (
        ch_input_bams,
        ch_cohort_varcalls
    )


    emit:
    ch_xhmm_varcalls = PER_SAMPLE_CALL_XHMM.out.xhmm_per_sample

}