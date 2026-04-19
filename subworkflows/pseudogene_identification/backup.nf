include { GETREFERENCE_FROM_REGION as GETREFERENCE_FROM_GENE_REGION} from "../../modules/pseudogene_identification/get_reference_from_region/main.nf"
include { GETREFERENCE_FROM_REGION as GETREFERENCE_FROM_PSEUDOGENE_REGION } from "../../modules/pseudogene_identification/get_reference_from_region/main.nf"
include { COMBINE_BAM } from "../../modules/pseudogene_identification/combine_bam/main.nf"
include { FORCEMAP } from "../../modules/pseudogene_identification/forcemap/main.nf"
include { DIFF_REFERENCE } from "../../modules/pseudogene_identification/diff_reference/main.nf"
include { READ_CLASSIFICATION } from "../../modules/pseudogene_identification/read_classification/main.nf"
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

    ch_gene_regions = ch_organism_gene_bed
        .splitCsv(header: false, sep: '\t')
        .combine(ch_reference_fasta)
        .map { row, reference -> tuple(reference, "${row[0]}:${row[1]}-${row[2]}") }


    GETREFERENCE_FROM_GENE_REGION(
        ch_gene_regions,
        ch_organism_gene_bed,
        ch_organism_pseudogene_bed
    )

    ch_gene_reference = GETREFERENCE_FROM_GENE_REGION.out[0]

    ch_pseudogene_regions = ch_organism_pseudogene_bed
        .splitCsv(header: false, sep: '\t')
        .combine(ch_reference_fasta)
        .map { row, reference -> tuple(reference, "${row[0]}:${row[1]}-${row[2]}") }

    GETREFERENCE_FROM_PSEUDOGENE_REGION(
        ch_pseudogene_regions,
        ch_organism_gene_bed,
        ch_organism_pseudogene_bed
    )



    ch_pseudogene_reference = GETREFERENCE_FROM_PSEUDOGENE_REGION.out[0]
    ch_combined_regions = ch_gene_reference.combine(ch_pseudogene_reference)
    ch_pseudogene_analysis = ch_input_bams
        .combine(ch_combined_regions)

    COMBINE_BAM(
        ch_pseudogene_analysis
    )

    ch_combined_bam = COMBINE_BAM.out.combined_bam

    FORCEMAP(
        ch_combined_bam
    )

    ch_diff_reference_input = ch_pseudogene_analysis
        .map { meta, bam, bai, gene_region, gene_ref, pseudogene_region, pseudogene_ref ->
            tuple(gene_region, pseudogene_region, gene_ref, pseudogene_ref)
        }

    DIFF_REFERENCE(
        ch_diff_reference_input
    )

    ch_read_processing_input = FORCEMAP.out.force_mapped_bam
        .combine(DIFF_REFERENCE.out.diff_output)
        .map { meta, gene_region, pseudogene_region, force_mapped_bam, dup_g_region, dup_pg_region, tsv ->
            tuple(meta, force_mapped_bam, gene_region, pseudogene_region, tsv)
        }
        .unique { meta, force_mapped_bam, gene_region, pseudogene_region, tsv ->
            "${meta.sample_name}|${gene_region}|${pseudogene_region}"
        }

    READ_PROCESSING(
        ch_read_processing_input
    )

    ch_read_classification_input = READ_PROCESSING.out.processed_reads

    READ_CLASSIFICATION(
        ch_read_classification_input
    )


}