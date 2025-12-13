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


    // Transform input channel
    ch_input_bams
        .map { sample ->
            sample[1..-1] // Remove the first element if necessary
        }
        .flatten() // Flatten the nested structure
        .collect() // Collect all files into a single list
        .set { ch_all_bam_files }

    // Debug transformed channel

    PREPARE_REGIONS_CNVKIT (
        ch_reference_fasta,
        ch_regions_of_interest,
        ch_all_bam_files
    )

    ch_prepared_regions = PREPARE_REGIONS_CNVKIT.out.prepared_regions

    GET_COVERAGE_CNVKIT (
        ch_input_bams,
        ch_prepared_regions
    )

    if (!params.normal_tumor) {
        ch_input_bams
            .count()
            .map { sample_count -> sample_count > 4 }
            .set { sample_number }
    } else {
        sample_number = Channel.value(true)
    }

    ch_coverage = GET_COVERAGE_CNVKIT.out.coverage

    // Collect all sample coverage files (target and antitarget from each sample)
    // Output is: tuple val(meta), path(target_cov), path(antitarget_cov)
    ch_sample_coverage_files = ch_coverage
        .map { meta, target_cov, antitarget_cov -> 
            [target_cov, antitarget_cov]  // Extract both coverage files
        }
        .collect()
        .map { pairs -> pairs.flatten() }  // Flatten after collecting to avoid duplicates

    // For tumor-only mode, pass empty list for normal coverage
    // For tumor-normal mode, the process will use the appropriate files
    ch_normal_coverage_files = Channel.value([])  // Empty list, not empty channel

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


    VARDICT_CALL (
        ch_input_bams,
        ch_reference_fasta,
        ch_reference_fasta_fai,
        ch_regions_of_interest
    )

    CNVKIT_CALL (
      ch_cnkvkit_segments
    )
    ch_cnvkit_varcall = CNVKIT_CALL.out.cnvkit_calls

    CNVKIT_CLASSIFY_AND_ANNOTATE(
        ch_cnvkit_varcall,
        ch_reference_gtf_tsv
    )

    emit:
    ch_vardict_vcfs = VARDICT_CALL.out.vcfs
}