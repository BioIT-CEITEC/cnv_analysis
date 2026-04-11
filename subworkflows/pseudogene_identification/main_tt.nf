include { EXTRACT_REGIONS } from "../../modules/pseudogene_identification/extract_regions/main.nf"
include { REALIGN_READS } from "../../modules/pseudogene_identification/realign_reads/main.nf"
include { READS_CLASSIFICATION } from "../../modules/pseudogene_identification/classify_reads/main.nf"


workflow PSEUDOGENE_ANALYSIS {

    take:
    ch_input_bams
    ch_reference_fasta
    ch_reference_fasta_fai
    ch_regions_bed

    main:

    EXTRACT_REGIONS(
        ch_input_bams,
        ch_regions_bed
    )

    ch_extracted_reads = EXTRACT_REGIONS.out.extracted_reads

    REALIGN_READS(
        ch_extracted_reads,
        ch_reference_fasta,
        ch_reference_fasta_fai
    )

    ch_realigned_reads = REALIGN_READS.out.realigned_reads

    READS_CLASSIFICATION(
        ch_realigned_reads,
        regions_bed
    )

}