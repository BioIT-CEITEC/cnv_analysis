include { PREPARE_REGIONS_CNVKIT } from "../../modules/cnvkit/prepare_regions/main.nf"
include { GET_COVERAGE_CNVKIT } from "../../modules/cnvkit/coverage/main.nf"
include { REFERENCE_CNVKIT } from "../../modules/cnvkit/reference_preparation/main.nf"
include { FIX_AND_SEGMENT_CNVKIT } from "../../modules/cnvkit/fix_and_segment/main.nf"
include { VARDICT_CALL } from "../../modules/vardict/main.nf"
include { CNVKIT_CALL } from "../../modules/cnvkit/call/main.nf"
include { DIAGRAM_AND_SCATTER_CNVKIT } from "../../modules/cnvkit/plotting/main.nf"

workflow CNVKIT_ANALYSIS {

    take:
    ch_input_bams
    ch_regions_of_interest
    ch_reference_fasta

    main:

    ch_input_bams
        .map { sample ->
            sample[1..-1]
        }
        .flatten()
        .set( ch_all_bam_files )

    PREPARE_REGIONS_CNVKIT (
        ch_reference_fasta
        ch_regions_of_interest
        ch_all_bam_files
    )

    ch_prepared_regions = PREPARE_REGIONS_CNVKIT.out.prepared_regions

    GET_COVERAGE_CNVKIT (
        ch_input_bams
        ch_prepared_regions
    )

    if (!params.tumor_normal) {
        ch_input_bams
            .size()
            .map { sample_count -> sample_count > 4 }
            .set { sample_number }
    } else {
        sample_number = true
    }

    ch_coverage = GET_COVERAGE_CNVKIT.out.coverage

    // Map tumor coverage files
    ch_tumor_coverage = ch_coverage.map { tuple -> 
        [tuple[1], tuple[2]] // Extract tumor target and antitarget files
    }

    // Map normal coverage files
    ch_normal_coverage = ch_coverage.map { tuple -> 
        [tuple[3], tuple[4]] // Extract normal target and antitarget files
    }

    // Flatten and collect tumor and normal coverage files into separate lists
    ch_tumor_coverage_files = ch_tumor_coverage.flatten().collect()
    ch_normal_coverage_files = ch_normal_coverage.flatten().collect()

    REFERENCE_CNVKIT (
        ch_reference_fasta,
        ch_tumor_coverage_files,
        ch_normal_coverage_files,
        ch_prepared_regions,
        sample_number
    )

    ch_cnvkit_ref = REFERENCE_CNVKIT.out.cnvkit_reference

    FIX_AND_SEGMENT_CNVKIT (
        ch_coverage
        ch_cnvkit_ref
    )

    ch_cnkvkit_segments = FIX_AND_SEGMENT_CNVKIT.out.cnvkit_segments

    VARDICT_CALL (
        ch_input_bams
        ch_reference_fasta
        ch_regions_of_interest
    )

    ch_vardict_vcfs = VARDICT_CALL.out.vcfs

    CNVKIT_CALL (
        ch_cnkvkit_segments
    )

    DIAGRAM_AND_SCATTER_CNVKIT (
        ch_input_bams
        ch_cnkvkit_segments
    )

}