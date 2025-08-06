include { AMBER } from "../../modules/amber/main.nf"
include { COBALT } from "../../modules/cobalt/main.nf"
include { PURPLE } from "../../modules/purple/main.nf"

workflow PURPLE_ANALYSIS {

  take:
  ch_input_bams
  ch_heterozygous_sites
  ch_region_bed
  ch_gc_profile
  ch_diploid_regions
  ch_vcfs
  ch_organism_fasta
  ch_organism_fasta_fai
  ch_organism_dict
  ch_organism_germline_hotspots
  ch_organism_driver_panel
  ch_organism_germline_dels
  ch_organism_vep
  
  main:
  
  AMBER (
    ch_input_bams,
    ch_heterozygous_sites,
    ch_region_bed
  )
  
  
  COBALT (
    ch_input_bams,
    ch_gc_profile,
    ch_diploid_regions,
    ch_region_bed
  )
  
  ch_purple_input = AMBER.out.amber_dir
    .join(COBALT.out.cobalt_dir)        // [ meta, amber, cobalt ]
    .join(ch_vcfs)           // [ meta, amber, cobalt, normal_vcf, tumor_vcf ]

  
  ch_purple_input.view()
  
  PURPLE (
    ch_purple_input,
    ch_organism_fasta,
    ch_organism_fasta_fai,
    ch_organism_dict,
    ch_gc_profile,
    ch_organism_germline_hotspots,
    ch_organism_driver_panel,
    ch_organism_vep,
    ch_organism_germline_dels
    
  )
  
}
