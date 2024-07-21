workflow CNVKIT_ANALYSIS {

    take:
    ch_input_bams
    ch_wgs_preprocess
    ch_reference_fasta

    main:

    PREPARE_REGIONS_CNVKIT (
        ch_reference_fasta
    )

    ch_input_bams
        .join(PREPARE_REGIONS_CNVKIT.out.cnvkit_beds, by: 'meta')
        .set { ch_input_cnvkit }

    GET_COVERAGE_CNVKIT (
        ch_input_cnvkit
    )

    ch_input_coverage = GET_COVERAGE_CNVKIT.out.coverage

    ch_normal_target = Channel.empty()
    ch_normal_antitarget = Channel.empty()

    if (params.normal_tumor) {
        ch_normal_target = ch_input_coverage.collect { it[1] }
        ch_normal_antitarget = ch_input_coverage.collect { it[2] }
}

    REFERENCE_CNVKIT (
        ch_reference_fasta
        ch_normal_target
        ch_normal_antitarget
    )

    ch_cnvkit_ref = REFERENCE_CNVKIT.out.reference

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