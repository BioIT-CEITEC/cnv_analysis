workflow JABCONTOOL_ANALYSIS {

    take:
    ch_input_bams
    ch_wgs_preprocess

    main:

    COVERAGE_CALC (
        ch_input_bams
        ch_wgs_preprocess
    )

    SNP_AF_CALC (
        ch_input_bams
    )

    SNP_AF_CALC.out
        .snpAF
        .collect { it[1] }
        .set { ch_input_tumor }

    if (params.normal_tumor) {
        SNP_AF_CALC.out
            .snpAF
            .collect { it[2] }
            .set { ch_input_normal }
    } else {
        ch_input_normal = Channel.empty()
    }

    JABCONTOOL_CALL (
        ch_input_tumor
        ch_input_normal
        ch_cohort_data
    )
}