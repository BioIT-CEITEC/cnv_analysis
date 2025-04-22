workflow CNVKIT_ANALYSIS {

    take:
    ch_input_bams
    ch_wgs_preprocess
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

def (ch_tumor_coverage, ch_normal_coverage) = GET_COVERAGE_CNVKIT.out.coverage
    .map { tuple ->
        def tumor = [tuple[1], tuple[2]]
        def normal = [tuple[3], tuple[4]]
        return [tumor, normal]
    }
    .transpose()

    ch_tumor_coverage_files = ch_tumor_coverage.flatten().collect()
    ch_normal_coverage_files = ch_normal_coverage.flatten().collect()

    REFERENCE_CNVKIT (
        ch_reference_fasta
        ch_tumor_coverage_files
        ch_normal_coverage_files
        ch_prepared_regions
        sample_number
    )

    ch_cnvkit_ref = REFERENCE_CNVKIT.out.cnvkit_reference

    VARDICT_CALL (
        ch_input_bams
    )

    ch_input_coverage
        .join(VARDICT_CALL.out.vcfs, by: 'meta')
        .set { ch_cnvkit_call_input }

    CNVKIT_CALL (
        ch_cnvkit_call_input
        ch_cnvkit_ref
    )
}