include { PREPARE_REGIONS_DELLY } from "../../modules/delly/prepare_regions/main.nf"
include { MERGE_REGIONS_DELLY } from "../../modules/delly/merge_regions/main.nf"
include { SV_CALLS_DELLY } from "../../modules/delly/sv_call/main.nf"
include { MERGE_AND_FILTER_SV_DELLY } from "../../modules/delly/merge_and_filter_sv/main.nf"
include { CNV_CALLS_DELLY } from "../../modules/delly/cnv_call/main.nf"

workflow DELLY_ANALYSIS {

    take:
    ch_input_bams
    ch_reference_fasta
    ch_reference_fasta_fai
    ch_excluded_regions
    ch_map_file

    main:
    
    PREPARE_REGIONS_DELLY(
    ch_input_bams,
    ch_reference_fasta,
    ch_reference_fasta_fai,
    ch_excluded_regions
    )
    
    ch_merge_regions_input = PREPARE_REGIONS_DELLY.out.prepared_regions_delly
      .map { meta, regions ->
        return regions }
      .collect()
      
    MERGE_REGIONS_DELLY(
    ch_merge_regions_input
    )
    
    ch_merged_regions = MERGE_REGIONS_DELLY.out.merged_regions_delly
    
    SV_CALLS_DELLY(
    ch_input_bams,
    ch_reference_fasta,
    ch_reference_fasta_fai,
    ch_merged_regions,
    ch_excluded_regions
    )
    
    ch_genotype_calls = SV_CALLS_DELLY.out.sv_genotype_delly
      .map { meta, calls, index -> 
        return calls
        }
      .collect()
    
    ch_genotype_index = SV_CALLS_DELLY.out.sv_genotype_delly
      .map { meta, calls, index -> 
        return index
        }
      .collect()
    
    MERGE_AND_FILTER_SV_DELLY(
    ch_genotype_calls,
    ch_genotype_index
    )
    
    ch_sv_calls = MERGE_AND_FILTER_SV_DELLY.out.final_sv_delly
    
    CNV_CALLS_DELLY(
    ch_input_bams,
    ch_reference_fasta,
    ch_reference_fasta_fai,
    ch_map_file,
    ch_sv_calls
    )
    
}