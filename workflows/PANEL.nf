 nextflow.enable.dsl = 2


// Validate inputs and potentially references

include { BED_PREPARATION } from "../modules/preprocessing/bed_preparation/main.nf"
include { BED_GTF_ANNOTATION } from "../modules/preprocessing/bed_gtf_annotation/main.nf"
include { PER_REGION_COVERAGE_CALC } from "../modules/preprocessing/coverage_preparation/main.nf"
include { COHORT_PREPARATION } from "../modules/preprocessing/cohort_preparation/main.nf"
include { CNVKIT_ANALYSIS } from "../subworkflows/cnvkit/main.nf"
include { JABCONTOOL_ANALYSIS } from "../subworkflows/jabcontool/main.nf"
include { GATK_ANALYSIS } from "../subworkflows/gatk/main.nf"
include { DELLY_ANALYSIS } from "../subworkflows/delly/main.nf"
include { CNMOPS_ANALYSIS } from "../subworkflows/cnMOPS/main.nf"
include { PANELCNMOPS_ANALYSIS } from "../subworkflows/panelcnMOPS/main.nf"
include { EXOMEDEPTH_ANALYSIS } from "../subworkflows/ExomeDepth/main.nf"
include { PSEUDOGENE_ANALYSIS } from "../subworkflows/pseudogene_identification/main.nf"
include { MERGE_VARIANT_CALLS } from "../modules/merge_and_smooth_CNVs_2/main.nf"
include { CLASSIFY_AND_ANNOTATE } from "../modules/classify_and_annotate_CNVs/main.nf"
include { ECOLE_ANALYSIS } from "../subworkflows/ECOLE/main.nf"
include { XHMM_ANALYSIS } from "../subworkflows/XHMM/main.nf"
include { CONIFER_ANALYSIS } from "../subworkflows/conifer/main.nf"
include { FREEC_ANALYSIS } from "../subworkflows/controlFREEC/main.nf"
include { CONSENSUS_MODEL_SCORE } from "../modules/consensus_model/main.nf"
include { CONSENSUS_MODEL_TRAIN } from "../modules/consensus_model/main.nf"
include { CONSENSUS_PREPARE_TRAINING_DIR } from "../modules/consensus_model/main.nf"
include { FINAL_FORMATTING_TABLES } from "../modules/final_formatting_tables/main.nf"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REFERENCE CONFIG LOADING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_organism_fasta = params.organism_fasta ? Channel.fromPath(params.organism_fasta).collect() : Channel.empty()
ch_organism_fasta_fai = params.organism_fasta ? Channel.fromPath(params.organism_fasta + '.fai').collect() : Channel.empty()
ch_organism_dict = params.organism_dict ? Channel.fromPath(params.organism_dict).collect() : Channel.empty()
ch_organism_dna_panel = params.organism_dna_panel ? Channel.fromPath(params.organism_dna_panel).collect() : Channel.empty()
ch_organism_cytoband = params.organism_cytoband ? Channel.fromPath(params.organism_cytoband).collect() : Channel.empty()
ch_organism_snps = params.organism_snps_panel ? Channel.fromPath(params.organism_snps_panel).collect() : Channel.empty()
ch_organism_ploidy_priors = params.organism_ploidy_priors ? Channel.fromPath(params.organism_ploidy_priors).collect() : Channel.empty()
ch_organism_gtf_tsv = params.organism_gtf_tsv ? Channel.fromPath(params.organism_gtf_tsv).collect() : Channel.empty()
ch_organism_gtf = params.organism_gtf ? Channel.fromPath(params.organism_gtf).collect() : Channel.empty()

Utils.loadSample(params)
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

    ch_gtf_for_coverage = params.organism_gtf
        ? ch_organism_gtf
        : Channel.value(file('/dev/null'))

    BED_GTF_ANNOTATION(
        ch_organism_dna_panel,
        ch_gtf_for_coverage
    )

    PER_REGION_COVERAGE_CALC(
        ch_samples,
        ch_organism_dna_panel,
        BED_GTF_ANNOTATION.out.annotation
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

    if (params.use_jabcontool) {
        ch_all_varcalls_for_merging = ch_all_varcalls
            .combine(JABCONTOOL_ANALYSIS.out.ch_jabcontool_norm_varcalls)
            .map { tuple ->
                def meta = tuple[0]
                def varcall_files = tuple[1]
                def jabcontool_file = tuple[2]
                def all_files = varcall_files + [jabcontool_file]
                [meta,[all_files].flatten()]
            }

    } else {
        ch_all_varcalls_for_merging = ch_all_varcalls
    }

    ch_coverage_for_final = PER_REGION_COVERAGE_CALC.out.region_coverage
        .map { meta, cov -> cov }
        .collect()

    // ─── Consensus model: optional training, then scoring ────────────────────
    // consensus_train = true trains a new model from ground-truth data before
    // scoring; otherwise the pre-trained model shipped with the module (or the
    // one given in params.consensus_model) is used.
    ch_consensus_model = Channel.empty()

    if (params.consensus_train) {
        if (!params.consensus_train_gt_dir) {
            error "consensus_train = true requires params.consensus_train_gt_dir (directory of <SAMPLE>_selected_regions.tsv ground-truth files)"
        }

        ch_consensus_gt_dir = Channel.fromPath(params.consensus_train_gt_dir, type: 'dir', checkIfExists: true).collect()

        // Training cohort: an external set of caller calls when
        // consensus_train_sv_dir is given, otherwise this run's own calls.
        if (params.consensus_train_sv_dir) {
            ch_consensus_train_dirs = Channel
                .fromPath("${params.consensus_train_sv_dir}/*", type: 'dir', checkIfExists: true)
                .collect()
        } else {
            CONSENSUS_PREPARE_TRAINING_DIR(ch_all_varcalls_for_merging)
            ch_consensus_train_dirs = CONSENSUS_PREPARE_TRAINING_DIR.out.sample_dir.collect()
        }

        CONSENSUS_MODEL_TRAIN(
            ch_consensus_train_dirs,
            ch_consensus_gt_dir
        )
        ch_consensus_model = CONSENSUS_MODEL_TRAIN.out.model.first()
    }
    else if (params.use_consensus) {
        ch_consensus_model = Channel.fromPath(
            params.consensus_model ?: "${projectDir}/modules/consensus_model/cnv_consensus_model.pkl",
            checkIfExists: true
        ).collect()
    }

    if (params.use_consensus) {
        CONSENSUS_MODEL_SCORE(
            ch_all_varcalls_for_merging,
            ch_consensus_model
        )
        ch_variants_to_annotate = CONSENSUS_MODEL_SCORE.out.merged_tsv
    } else {
        MERGE_VARIANT_CALLS(
            ch_all_varcalls_for_merging,
            ch_organism_dna_panel,
            ch_organism_gtf
        )
        ch_variants_to_annotate = MERGE_VARIANT_CALLS.out.merged_tsv
    }

    CLASSIFY_AND_ANNOTATE(
        ch_variants_to_annotate,
        ch_organism_gtf_tsv
    )

    ch_for_final = ch_variants_to_annotate
        .map { meta, tsv, bed -> tsv }
        .mix(
            CLASSIFY_AND_ANNOTATE.out.annotated_tsv
                .map { meta, tsv -> tsv }
        )
        .collect()

    FINAL_FORMATTING_TABLES(
        ch_for_final,
        ch_organism_dna_panel,
        ch_organism_gtf,
        ch_coverage_for_final
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