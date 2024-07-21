include { AMBER_BAF_CALCULATION } from './modules/amber_baf_calculation.nf'
include { COBALT_READ_DEPTH } from './modules/cobalt_read_depth.nf'
include { GRIDSS_CALL } from './modules/gridss_call.nf'
include { GRIPSS_CALL } from './modules/gripss_call.nf'
include { PURPLE_CALL } from './modules/purple_call.nf'

workflow PURPLE_ANALYSIS {

    take:
    ch_input_bams

    main:

    AMBER_BAF_CALCULATION {
        ch_input_bams
    }

    COBALT_READ_DEPTH {
        ch_input_bams
    }

    GRIDSS_CALL {
        ch_input_bams
    }

    ch_input_vcf = GRIDSS_CALL.out.vcfs

    GRIPSS_CALL {
        ch_input_vcf
    }

    AMBER_BAF_CALCULATION.out.
        .amber
        .join(COBALT_READ_DEPTH.out.cobalt, by: 'meta')
        .join(GRIPSS_CALL.out.vcfs, by: 'meta')
        .set { ch_input_purple }

    PURPLE_CALL {
        ch_input_purple
    }

    emit:
    ch_purple_outputs = PURPLE_CALL.out.cnvs
}