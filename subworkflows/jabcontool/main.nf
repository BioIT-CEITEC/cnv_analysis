include { COVERAGE_CALC } from "../../modules/jabcontool/coverage/main.nf"
include { SNP_AF_CALC } from "../../modules/jabcontool/AF_calculation/main.nf"
include { JABCONTOOL_CALL } from "../../modules/jabcontool/call/main.nf"
include { VARIANT_NORMALIZATION_JABCONTOOL } from "../../modules/jabcontool/variantNormalization/main.nf"

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

    COVERAGE_CALC (
        ch_input_bams,
        ch_organism_dna_panel,
        ch_reference_fasta_fai
    )
    SNP_AF_CALC(
        ch_input_bams,
        ch_reference_fasta,
        ch_reference_fasta_fai,
        ch_organism_snps
    )
    ch_combined = COVERAGE_CALC.out.region_coverage
        .map { meta, cov -> cov }
        .collect()
        .map { cov_files -> [cov_files] }
        .concat(
            SNP_AF_CALC.out.snpAF
                .map { meta, snp -> snp }
                .collect()
                .map { snp_files -> [snp_files] }
        )
        .collect()
        .map { lists -> tuple(lists[0], lists[1]) }

    JABCONTOOL_CALL (
        ch_combined,
        ch_organism_dna_panel,
        ch_gc_profile,
        ch_organism_cytoband,
        ch_organism_snps
        )

    ch_jabcontool_varcalls = JABCONTOOL_CALL.out.final_CNV_probs

    VARIANT_NORMALIZATION_JABCONTOOL (
        ch_jabcontool_varcalls
    )
    emit:
    ch_jabcontool_norm_varcalls = VARIANT_NORMALIZATION_JABCONTOOL.out.normalized_varcalls
}
