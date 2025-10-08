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
    COVERAGE_CALC.out.region_coverage
    .join( SNP_AF_CALC.out.snpAF )            // single synchronisation
    .multiMap { meta, normal_cov, normal_snp ->

        // one key/value pair per output channel
 
        normal_cov  : normal_cov
        normal_snps : normal_snp
    }
    .set { ch_all }  


def normal_cov_collected = ch_all.normal_cov.toList()
def normal_snps_collected = ch_all.normal_snps.toList()

def combined_channel = normal_cov_collected
    .concat(normal_snps_collected)
    .concat(normal_snps_collected)
    .toList()
    .map { normal_covs, normal_snps ->
        return [normal_covs, normal_snps]
    }

    JABCONTOOL_CALL (
        combined_channel,
        ch_region_bed,
        ch_gc_profile,
        ch_organism_cytoband,
        ch_organism_snps
        )
}