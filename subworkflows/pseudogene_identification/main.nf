include { ALIGN_REGIONS } from "../../modules/pseudogene_identification/align_regions/main.nf"
include { EXTRACT_REGIONS } from "../../modules/pseudogene_identification/extract_regions/main.nf"
include { REALIGN_READS } from "../../modules/pseudogene_identification/realign_reads/main.nf"
include { CLASSIFY_READS as READS_CLASSIFICATION } from "../../modules/pseudogene_identification/classify_reads/main.nf"
include { MISMATCH_PILEUP } from "../../modules/pseudogene_identification/mismatch_pileup/main.nf"

workflow PSEUDOGENE_ANALYSIS {

    take:
    ch_input_bams
    ch_reference_fasta
    ch_reference_fasta_fai
    ch_regions_bed

    main:

    ALIGN_REGIONS(
        ch_reference_fasta,
        ch_reference_fasta_fai,
        ch_regions_bed
    )

    EXTRACT_REGIONS(
        ch_input_bams,
        ch_regions_bed
    )

    ch_extracted_reads = EXTRACT_REGIONS.out.extracted_reads

    REALIGN_READS(
        ch_extracted_reads,
        ch_reference_fasta,
        ch_reference_fasta_fai,
        ch_regions_bed
    )

    ch_realigned_reads = REALIGN_READS.out.realigned_reads

    MISMATCH_PILEUP(
        ch_realigned_reads,
        ALIGN_REGIONS.out.diff_tsvs.first(),
        ch_regions_bed
    )

    READS_CLASSIFICATION(
        ch_realigned_reads,
        ALIGN_REGIONS.out.diff_tsvs.first(),
        ch_regions_bed
    )

    emit:
    pileup_tsvs      = MISMATCH_PILEUP.out.pileup_tsv
    classified_reads = READS_CLASSIFICATION.out.classified_reads

}