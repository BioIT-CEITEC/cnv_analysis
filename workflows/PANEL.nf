nextflow.enable.dsl = 2


// Validate inputs and potentially references

include { BED_PREPARATION } from "../modules/preprocessing/bed_preparation/main.nf"
include { COHORT_PREPARATION } from "../modules/preprocessing/cohort_preparation/main.nf"
include { CNVKIT_ANALYSIS } from "../subworkflows/cnvkit/main.nf"
include { JABCONTOOL_ANALYSIS } from "../subworkflows/jabcontool/main.nf"
include { GATK_ANALYSIS } from "../subworkflows/gatk/main.nf"
include { DELLY_ANALYSIS } from "../subworkflows/delly/main.nf"
include { CNMOPS_ANALYSIS } from "../subworkflows/cnMOPS/main.nf"
include { PANELCNMOPS_ANALYSIS } from "../subworkflows/panelcnMOPS/main.nf"
include { EXOMEDEPTH_ANALYSIS } from "../subworkflows/ExomeDepth/main.nf"
include { VARIANT_NORMALIZATION } from "../modules/variantNormalization/main.nf"
include { PSEUDOGENE_ANALYSIS } from "../subworkflows/pseudogene_identification/main.nf"
include { MERGE_VARIANT_CALLS } from "../modules/merge_and_smooth_CNVs/main.nf"
include { CLASSIFY_AND_ANNOTATE } from "../modules/classify_and_annotate_CNVs/main.nf"
include { ECOLE_ANALYSIS } from "../subworkflows/ECOLE/main.nf"
include { XHMM_ANALYSIS } from "../subworkflows/XHMM/main.nf"
include { CONIFER_ANALYSIS } from "../subworkflows/conifer/main.nf"
include { FREEC_ANALYSIS } from "../subworkflows/controlFREEC/main.nf"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REFERENCE CONFIG LOADING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Load the reference data from the paths in config file
params = Utils.load_organism(params)
//Utils.set_data_tags(params)
Utils.loadSample(params)

ch_organism_fasta = params.organism_fasta ? Channel.fromPath(params.organism_fasta).collect() : Channel.empty()
ch_organism_fasta_fai = params.organism_fasta ? Channel.fromPath(params.organism_fasta + '.fai').collect() : Channel.empty()
ch_organism_dict = params.organism_dict ? Channel.fromPath(params.organism_dict).collect() : Channel.empty()
ch_organism_dna_panel = params.organism_dna_panel ? Channel.fromPath(params.organism_dna_panel).collect() : Channel.empty()
ch_organism_cytoband = params.organism_cytoband ? Channel.fromPath(params.organism_cytoband).collect() : Channel.empty()
ch_organism_snps = params.organism_snps_panel ? Channel.fromPath(params.organism_snps_panel).collect() : Channel.empty()
ch_heterozygous_sites = params.organism_hetsites ? Channel.fromPath(params.organism_hetsites).collect() : Channel.empty()
ch_gc_content = params.organism_gc_profile ? Channel.fromPath(params.organism_gc_profile).collect() : Channel.empty()
ch_diploid_regions = params.organism_diploid_regions ? Channel.fromPath(params.organism_diploid_regions).collect() : Channel.empty()
ch_organism_germline_hotspots = params.organism_germline_hotspots ? Channel.fromPath(params.organism_germline_hotspots).collect() : Channel.empty()
ch_organism_driver_panel = params.organism_germline_driverpanel ? Channel.fromPath(params.organism_germline_driverpanel).collect() : Channel.empty()
ch_organism_germline_dels = params.organism_germline_dels ? Channel.fromPath(params.organism_germline_dels).collect() : Channel.empty()
ch_organism_vep = params.organism_vep_dir ? Channel.fromPath(params.organism_vep_dir).collect() : Channel.empty()
ch_organism_ploidy_priors = params.organism_ploidy_priors ? Channel.fromPath(params.organism_ploidy_priors).collect() : Channel.empty()
ch_organism_excluded_sites = params.organism_excluded_sites ? Channel.fromPath(params.organism_excluded_sites).collect() : Channel.empty()
ch_organism_delly_map = params.organism_delly_map ? Channel.fromPath(params.organism_delly_map).collect() : Channel.empty()
ch_organism_gtf_tsv = params.organism_gtf_tsv ? Channel.fromPath(params.organism_gtf_tsv).collect() : Channel.empty()
ch_organism_gene_bed = params.organism_gene_bed ? Channel.fromPath(params.organism_gene_bed).collect() : Channel.empty()
ch_organism_pseudogene_bed = params.organism_pseudogene_bed ? Channel.fromPath(params.organism_pseudogene_bed).collect() : Channel.empty()
ch_organism_gtf = params.organism_gtf ? Channel.fromPath(params.organism_gtf).collect() : Channel.empty()

def inputData = Utils.parseInputVC(params.new_samples, projectDir, log)

def samples = inputData.samples

workflow PANEL_WES {

   ch_cohort_data = Channel.empty()
   ch_all_varcalls = Channel.empty()

  ch_samples = Channel.fromList(samples)
    .map { meta, bamPath, baiPath -> tuple(meta, file(bamPath), file(baiPath)) }

    BED_PREPARATION (
        ch_organism_fasta,
        ch_organism_dict
    )

    ch_binned_genome = BED_PREPARATION.out.binned_genome
    ch_gc_profile = BED_PREPARATION.out.gc_profile

    if (params.use_cnvkit) {
        CNVKIT_ANALYSIS (
        ch_samples,
        ch_organism_dna_panel,
        ch_organism_fasta,
        ch_organism_fasta_fai,
        ch_organism_gtf_tsv
        )

        ch_all_varcalls = Utils.mixAndCollectVarcalls(
            ch_all_varcalls,
            CNVKIT_ANALYSIS.out.ch_cnvkit_calls.map { tuple -> [tuple[0], tuple[5]] }
        )
    }

    if (params.use_jabcontool) {
        JABCONTOOL_ANALYSIS (
        ch_samples,
        ch_organism_fasta,
        ch_organism_fasta_fai,
        ch_organism_snps,
        ch_organism_dna_panel,
        ch_organism_cytoband,
        ch_gc_profile,
        ch_cohort_data
        )
    }

    if (params.use_gatk) {
        GATK_ANALYSIS(
        ch_samples,
        ch_organism_fasta,
        ch_organism_fasta_fai,
        ch_organism_dna_panel,
        ch_organism_dict,
        ch_organism_ploidy_priors
        )

        ch_all_varcalls = Utils.mixAndCollectVarcalls(
            ch_all_varcalls,
            GATK_ANALYSIS.out.ch_gatk_segment.map { tuple -> [tuple[0], tuple[1]] }
        )
}

    if (params.use_ecole) {
        ECOLE_ANALYSIS(
            ch_samples,
            ch_organism_dna_panel
        )
        ch_all_varcalls = Utils.mixAndCollectVarcalls(ch_all_varcalls, ECOLE_ANALYSIS.out.ch_ecole_varcalls)
    }

    if (params.use_xhmm) {
        XHMM_ANALYSIS(
            ch_samples,
            ch_organism_dna_panel
        )
        ch_all_varcalls = Utils.mixAndCollectVarcalls(ch_all_varcalls, XHMM_ANALYSIS.out.ch_xhmm_varcalls)
    }

    if (params.use_conifer) {
        CONIFER_ANALYSIS(
            ch_samples,
            ch_organism_dna_panel,
            ch_organism_gtf
        )
        ch_all_varcalls = Utils.mixAndCollectVarcalls(ch_all_varcalls, CONIFER_ANALYSIS.out.ch_conifer_varcalls)
    }

    if (params.use_controlfreec) {
        FREEC_ANALYSIS(
            ch_samples,
            ch_organism_fasta,
            ch_organism_fasta_fai,
            ch_organism_dna_panel
        )
        ch_all_varcalls = Utils.mixAndCollectVarcalls(ch_all_varcalls, FREEC_ANALYSIS.out.ch_freec_varcalls)
    }

    ch_samples_as_cohort = ch_samples
        .map { meta, bam, bai -> bam }
        .collect()
        .map { bam_files -> [bam_files] }
        .concat(
            ch_samples
                .map { meta, bam, bai -> bai }
                .collect()
                .map { bai_files -> [bai_files] }
        )
        .collect()
        .map { lists -> tuple(lists[0], lists[1]) }

    ch_panelcnmops_cohort = params.use_panelcnmops_cohortdata
        ? Channel.fromPath(params.cohort_panelcnmops_ref).collect()
        : ch_samples_as_cohort

  if (params.use_panelcnmops) {
        PANELCNMOPS_ANALYSIS (
        ch_samples,
        ch_panelcnmops_cohort,
        ch_organism_dna_panel
        )

        ch_all_varcalls = Utils.mixAndCollectVarcalls(ch_all_varcalls, PANELCNMOPS_ANALYSIS.out.ch_panelcnmops_varcalls)
  }

    ch_cnmops_cohort = params.use_cnmops_cohortdata
        ? Channel.fromPath(params.cohort_cnmops_ref).collect()
        : ch_samples_as_cohort

  if (params.use_cnmops) {
        CNMOPS_ANALYSIS (
        ch_samples,
        ch_cnmops_cohort,
        ch_organism_dna_panel
        )

        ch_all_varcalls = Utils.mixAndCollectVarcalls(ch_all_varcalls, CNMOPS_ANALYSIS.out.ch_cnmops_varcalls)

  }
    ch_exomedepth_cohort = params.use_exomedepth_cohortdata
        ? Channel.fromPath(params.cohort_exomedepth_ref).collect()
        : ch_samples_as_cohort

  if (params.use_exomedepth) {
        EXOMEDEPTH_ANALYSIS (
        ch_samples,
        ch_exomedepth_cohort,
        ch_organism_dna_panel,
        ch_organism_fasta
        )

        ch_all_varcalls = Utils.mixAndCollectVarcalls(ch_all_varcalls, EXOMEDEPTH_ANALYSIS.out.ch_exomeDepth_varcalls)

  }

    VARIANT_NORMALIZATION(
        ch_all_varcalls
    )

    if (params.use_jabcontool) {
        ch_all_varcalls_for_merging = VARIANT_NORMALIZATION.out.normalized_varcalls
            .combine(JABCONTOOL_ANALYSIS.out.ch_jabcontool_norm_varcalls)
            .map { tuple ->
                def meta = tuple[0]
                def varcall_files = tuple[1]
                def jabcontool_file = tuple[2]
                def all_files = varcall_files + [jabcontool_file]
                [meta,[all_files].flatten()]
            }

    } else {
        ch_all_varcalls_for_merging = VARIANT_NORMALIZATION.out.normalized_varcalls
    }

    MERGE_VARIANT_CALLS(
        ch_all_varcalls_for_merging,
        ch_organism_dna_panel
    )

    ch_variants_to_annotate = MERGE_VARIANT_CALLS.out.merged_calls

    CLASSIFY_AND_ANNOTATE(
        ch_variants_to_annotate,
        ch_organism_gtf_tsv
    )

    if (params.run_pseudogene) {
        PSEUDOGENE_ANALYSIS(
            ch_samples,
            ch_organism_fasta,
            ch_organism_fasta_fai,
            ch_organism_gene_bed
        )
    }
}