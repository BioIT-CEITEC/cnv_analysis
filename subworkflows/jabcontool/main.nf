include { COVERAGE_CALC } from "../../modules/jabcontool/coverage/main.nf"
include { SNP_AF_CALC } from "../../modules/jabcontool/AF_calculation/main.nf"
include { JABCONTOOL_CALL } from "../../modules/jabcontool/call/main.nf"

workflow JABCONTOOL_ANALYSIS {

    take:
    ch_input_bams
    ch_reference_fasta
    ch_reference_fasta_fai
    ch_organism_snps
    ch_organism_dna_panel
    ch_organism_cytoband
    ch_cohort_data

    main:
    ch_region_bed = Channel.of(params.tumor_normal_paired ? ch_binned_genome : ch_organism_dna_panel)

    COVERAGE_CALC (
        ch_input_bams
        ch_region_bed
        ch_reference_fasta_fai
    )

    ch_coverage = COVERAGE_CALC.out.region_coverage

    SNP_AF_CALC(
        ch_input_bams
        ch_reference_fasta
        ch_reference_fasta_fai
        ch_organism_snps
    )

     ch_snps = SNP_AF_CALC.out.snp_coverage

    ch_cov_snps = ch_coverage
        .join(ch_snps, by: 'meta')
        .map { tuple_coverage, tuple_snps ->
            def tumor_cov = tuple_coverage[1]
            def normal_cov = tuple_coverage[2]
            def tumor_snps = tuple_snps[1]
            def normal_snps = tuple_snps[2]
            return [tumor_cov, tumor_snps, normal_cov, normal_snps]
        }
        .groupTuple()

    JABCONTOOL_CALL(
        ch_cov_snps
        ch_region_bed
        ch_gc_profile
        ch_organism_cytoband
        ch_cohort_data
        ch_organism_snps
    )






}