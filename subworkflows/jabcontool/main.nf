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
    ch_gc_profile
    ch_cohort_data

    main:
    def tumorNormal = params.tumor_normal
    ch_region_bed = tumorNormal ? ch_binned_genome : ch_organism_dna_panel

    COVERAGE_CALC (
        ch_input_bams,
        ch_region_bed,
        ch_reference_fasta_fai
    )
    SNP_AF_CALC(
        ch_input_bams,
        ch_reference_fasta,
        ch_reference_fasta_fai,
        ch_organism_snps
    )
    // Collect all coverage and SNP files and create a proper tuple structure
    ch_combined = COVERAGE_CALC.out.region_coverage
        .map { meta, cov -> cov }
        .collect()
        .map { cov_files -> [cov_files] }  // Wrap in list to preserve structure
        .concat(
            SNP_AF_CALC.out.snpAF
                .map { meta, snp -> snp }
                .collect()
                .map { snp_files -> [snp_files] }  // Wrap in list to preserve structure
        )
        .collect()
        .map { lists -> tuple(lists[0], lists[1]) }  // Create tuple from the two lists

    JABCONTOOL_CALL (
        ch_combined,
        ch_region_bed,
        ch_gc_profile,
        ch_organism_cytoband,
        ch_organism_snps
        )
}