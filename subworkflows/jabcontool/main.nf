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
    ch_region_bed = params.normal_tumor ? ch_binned_genome : ch_organism_dna_panel

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
    .multiMap { meta, tumor_cov, normal_cov, tumor_snp, normal_snp ->

        // one key/value pair per output channel
        tumor_cov   : tumor_cov 
        tumor_snps  : tumor_snp
        normal_cov  : normal_cov
        normal_snps : normal_snp
    }
    .set { ch_all }  

        
def tumor_cov_collected = ch_all.tumor_cov.toList()
def tumor_snps_collected = ch_all.tumor_snps.toList()
def normal_cov_collected = ch_all.normal_cov.toList() 
def normal_snps_collected = ch_all.normal_snps.toList()

def combined_channel = tumor_cov_collected
    .concat(tumor_snps_collected)
    .concat(normal_cov_collected)
    .concat(normal_snps_collected)
    .toList()
    .map { tumor_covs, tumor_snps, normal_covs, normal_snps ->
        return [tumor_covs, tumor_snps, normal_covs, normal_snps]
    }
    
    JABCONTOOL_CALL (
        combined_channel,
        ch_region_bed,
        ch_gc_profile,
        ch_organism_cytoband,
        ch_organism_snps
        
        )
    
}