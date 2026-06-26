include { GET_PROBES_CONIFER } from "../../modules/conifer/get_probes/main.nf"
include { CALCULATE_RPKM_CONIFER } from "../../modules/conifer/calculate_rpkm/main.nf"
include { COHORT_ANALYSIS_CONIFER } from "../../modules/conifer/cohort_analysis/main.nf"
include { COHORT_CALL_CONIFER } from "../../modules/conifer/cohort_call/main.nf"
include { EXPORT_SAMPLE_CONIFER } from "../../modules/conifer/export_sample/main.nf"
include { EXTRACT_SAMPLE_CONIFER } from "../../modules/conifer/extract_sample/main.nf"

workflow CONIFER_ANALYSIS {

    take:
    ch_input_bams
    ch_regions_of_interest
    ch_gtf_file

    main:

    ch_input_bams
        .map { sample ->
            sample[1..-1]
        }
        .flatten()
        .collect()
        .set { ch_all_bam_files }

    GET_PROBES_CONIFER(
        ch_regions_of_interest,
        ch_gtf_file
    )

    ch_region_probes = GET_PROBES_CONIFER.out.conifer_probes

    CALCULATE_RPKM_CONIFER (
        ch_input_bams,
        ch_region_probes
    )

    CALCULATE_RPKM_CONIFER.out.rpkm_conifer
        .map { sample, rpkm -> rpkm }
        .collect()
        .set { ch_all_rpkm_files }

    COHORT_ANALYSIS_CONIFER (
        ch_all_rpkm_files,
        ch_region_probes
    )

    ch_conifer_cohort_analysis = COHORT_ANALYSIS_CONIFER.out.conifer_cohort_analysis

    COHORT_CALL_CONIFER (
        ch_conifer_cohort_analysis
    )

    ch_conifer_cohort_calls = COHORT_CALL_CONIFER.out.conifer_cohort_calls

    ch_input_bams
    .map { meta, bam, bai -> meta }
    .combine(ch_conifer_cohort_analysis)
    .set { ch_meta_hdf5 }

    EXPORT_SAMPLE_CONIFER (
        ch_meta_hdf5
    )

    EXTRACT_SAMPLE_CONIFER (
        ch_input_bams,
        ch_conifer_cohort_calls
    )

    emit:
    ch_conifer_varcalls = EXTRACT_SAMPLE_CONIFER.out.conifer_per_sample

}