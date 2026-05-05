include { PREPARE_REGIONS_CNVKIT } from "../../modules/cnvkit/prepare_regions/main.nf"
include { GET_COVERAGE_CNVKIT } from "../../modules/cnvkit/coverage/main.nf"
include { REFERENCE_CNVKIT } from "../../modules/cnvkit/reference_preparation/main.nf"
include { FIX_AND_SEGMENT_CNVKIT } from "../../modules/cnvkit/fix_and_segment/main.nf"
include { VARDICT_CALL } from "../../modules/vardict/main.nf"
include { CNVKIT_CALL } from "../../modules/cnvkit/call/main.nf"
include { DIAGRAM_AND_SCATTER_CNVKIT } from "../../modules/cnvkit/plotting/main.nf"
include { CNVKIT_CLASSIFY_AND_ANNOTATE } from "../../modules/cnvkit/classify_and_annotate/main.nf"

workflow CNVKIT_ANALYSIS {

    take:
    ch_input_bams
    ch_regions_of_interest
    ch_reference_fasta
    ch_reference_fasta_fai
    ch_reference_gtf_tsv

    main:

    PREPARE_REGIONS_CNVKIT (
        ch_reference_fasta,
        ch_regions_of_interest
    )

    ch_prepared_regions = PREPARE_REGIONS_CNVKIT.out.prepared_regions

    GET_COVERAGE_CNVKIT (
        ch_input_bams,
        ch_prepared_regions
    )

    ch_input_bams
        .count()
        .map { sample_count -> sample_count > 4 }
        .set { sample_number }

    ch_coverage = GET_COVERAGE_CNVKIT.out.coverage

    ch_sample_coverage_files = ch_coverage
        .map { meta, target_cov, antitarget_cov -> 
            [target_cov, antitarget_cov]
        }
        .collect()
        .map { pairs -> pairs.flatten() }

    ch_normal_coverage_files = Channel.value([])

    REFERENCE_CNVKIT (
        ch_reference_fasta,
        ch_sample_coverage_files,
        ch_normal_coverage_files,
        ch_prepared_regions,
        sample_number
    )

    ch_cnvkit_ref = REFERENCE_CNVKIT.out.cnvkit_reference

    FIX_AND_SEGMENT_CNVKIT(
        ch_coverage,
        ch_cnvkit_ref
    )

    ch_cnkvkit_segments = FIX_AND_SEGMENT_CNVKIT.out.cnvkit_segments

    CNVKIT_CALL (
      ch_cnkvkit_segments
    )

    emit:
    ch_cnvkit_calls = CNVKIT_CALL.out.cnvkit_calls
}