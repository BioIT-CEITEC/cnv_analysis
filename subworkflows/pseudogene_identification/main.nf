include { GETREFERENCE_FROM_REGION } from "../../modules/pseudogene_identification/get_reference_from_region/main.nf"
include { COMBINE_BAM } from "../../modules/pseudogene_identification/combine_bam/main.nf"
include { FORCEMAP } from "../../modules/pseudogene_identification/forcemap/main.nf"
include { DIFF_REFERENCE } from "../../modules/pseudogene_identification/diff_reference/main.nf"
include { READ_CLASSIFICTION } from "../../modules/pseudogene_identification/read_classification/main.nf"
include { READ_PROCESSING } from "../../modules/pseudogene_identification/read_processing/main.nf"
include { SPLIT_BAM } from "../../modules/pseudogene_identification/split_bam/main.nf"
include { SUMMARIZE } from "../../modules/pseudogene_identification/summarize/main.nf"

workflow PSEUDOGENE_ANALYSIS {

    take:
    ch_input_bams
    ch_reference_fasta
    ch_reference_fasta_fai
    ch_organism_gene_bed
    ch_organism_pseudogene_bed

    main:

    gene_regions_ch = ch_organism_gene_bed
        .splitCsv(header: false, sep: '\t')
        .combine(ch_reference_fasta)
        .map { row, reference -> tuple(reference, "${row[0]}:${row[1]}-${row[2]}") }

    gene_regions_ch.view()

    //GETREFERENCE_FROM_REGION(
    //    gene_regions_ch
    //)

}
