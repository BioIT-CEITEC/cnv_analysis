include { GENMAP_INDEX_FREEC } from "../../modules/control_freec/genmap_index/main.nf"
include { GENMAP_MAP_FREEC } from "../../modules/control_freec/genmap_map/main.nf"
include { CNV_CALL_FREEC } from "../../modules/control_freec/cnv_call/main.nf"
include { EXTRACT_SAMPLE_FREEC } from "../../modules/control_freec/extraction_sample/main.nf"

workflow FREEC_ANALYSIS {

    take:
    ch_input_bams
    ch_reference_fasta
    ch_reference_fasta_fai
    ch_bed_file 


    main:

    GENMAP_INDEX_FREEC (
        ch_reference_fasta,
        ch_reference_fasta_fai
    )

    ch_genmap_index = GENMAP_INDEX_FREEC.out.genmap_index
    ch_freec_chrlen = GENMAP_INDEX_FREEC.out.freec_chrlen

    GENMAP_MAP_FREEC (
        ch_genmap_index
    )

    ch_mappability_bg = GENMAP_MAP_FREEC.out.mappability_bg

    CNV_CALL_FREEC (
        ch_input_bams,
        ch_freec_chrlen,
        ch_mappability_bg,
        ch_bed_file
    )

    ch_freec_calls = CNV_CALL_FREEC.out.var_call

     EXTRACT_SAMPLE_FREEC (
        ch_freec_calls
    )

    emit:
    ch_freec_varcalls = EXTRACT_SAMPLE_FREEC.out.cnv_tsv

}