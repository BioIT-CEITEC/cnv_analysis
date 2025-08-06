include { PREPARE_REGIONS_GATK } from "../../modules/gatk/prepare_regions/main.nf"
include { ANNOTATE_REGIONS_GATK } from "../../modules/gatk/annotate_regions/main.nf"
include { COUNT_READS_GATK } from "../../modules/gatk/count_reads/main.nf"
include { FILTER_INTERVALS_GATK } from "../../modules/gatk/filter_intervals/main.nf"
include { GERMLINE_PLOIDY_DETERMINATION_GATK } from "../../modules/gatk/ploidy_determination/main.nf"
include { GERMLINE_VARCALLS_GATK } from "../../modules/gatk/cnv_varcalls/main.nf"
include { POSTPROCESSING_CNV_GATK } from "../../modules/gatk/cnv_postprocessing/main.nf"

workflow GATK_ANALYSIS {

    take:
    normal_bam
    organism_fasta
    organism_fasta_fai
    organism_regions
    organism_dict
    organism_ploidy_priors
    
    main:
    
    PREPARE_REGIONS_GATK (
    organism_fasta,
    organism_fasta_fai,
    organism_regions,
    organism_dict
    )
    
    ch_prepared_regions = PREPARE_REGIONS_GATK.out.prepared_regions_gatk.collect()
    
    ANNOTATE_REGIONS_GATK (
    organism_fasta,
    organism_fasta_fai,
    organism_regions,
    organism_dict,
    ch_prepared_regions
    )
    
    ch_annotated_regions = ANNOTATE_REGIONS_GATK.out.annotated_regions_gatk.collect()
    ch_annotated_regions.view()
    
    
    COUNT_READS_GATK(
    normal_bam,
    organism_fasta,
    organism_fasta_fai,
    organism_regions,
    organism_dict,
    ch_prepared_regions
    )
   
   
    ch_filtering_input = COUNT_READS_GATK.out.read_counts_gatk
      .map { meta, counts -> 
        return counts}
      .collect()
      
    FILTER_INTERVALS_GATK(
    ch_prepared_regions,
    ch_annotated_regions,
    organism_regions,
    ch_filtering_input
    )
    
    ch_qc_filtered_intervals = FILTER_INTERVALS_GATK.out.qc_filtered_regions_gatk
    
    GERMLINE_PLOIDY_DETERMINATION_GATK(
    ch_qc_filtered_intervals,
    ch_filtering_input,
    organism_ploidy_priors
    )
    
    ch_ploidy_determination = GERMLINE_PLOIDY_DETERMINATION_GATK.out.ploidy_model_gatk.collect()
    
    GERMLINE_VARCALLS_GATK(
    ch_prepared_regions,
    ch_filtering_input,
    ch_annotated_regions,
    ch_ploidy_determination
    )
    
    ch_cohort_data = GERMLINE_VARCALLS_GATK.out.germline_varcalls_gatk.collect()
  
      
    POSTPROCESSING_CNV_GATK(
    normal_bam,
    ch_cohort_data,
    ch_ploidy_determination,
    organism_dict
    )
    
    ch_gatk_segment = POSTPROCESSING_CNV_GATK.out.postprocess_cnv_gatk
      .map { meta, intervals_vcf, intervals_index, segments_vcf, segments_index, copy_ratios -> 
        return[meta, segments_vcf]
      }
    
    emit:
    ch_gatk_segment

}