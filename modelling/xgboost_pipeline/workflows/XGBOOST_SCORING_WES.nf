nextflow.enable.dsl = 2

include { XGBOOST_CONSENSUS_SCORE } from '../modules/xgboost_consensus_model/main.nf'
include { CLASSIFY_AND_ANNOTATE }   from '../../../modules/classify_and_annotate_CNVs/main.nf'
include { FINAL_FORMATTING_TABLES } from '../../../modules/final_formatting_tables/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REFERENCE CONFIG LOADING (same params.organism_* as the main pipeline)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_organism_dna_panel = params.organism_dna_panel ? Channel.fromPath(params.organism_dna_panel).collect() : Channel.empty()
ch_organism_gtf       = params.organism_gtf ? Channel.fromPath(params.organism_gtf).collect() : Channel.empty()
ch_organism_gtf_tsv   = params.organism_gtf_tsv ? Channel.fromPath(params.organism_gtf_tsv).collect() : Channel.empty()

// Same meta/path-defaulting convention as Utils.parseInputVC's BAM-path
// defaulting (lib/Utils.groovy), adapted for pre-computed caller-output
// directories instead of BAMs.
def parseXgbSamples(samplesMap, projectDir) {
    return samplesMap.collect { key, value ->
        def meta = [sample_name: value.sample_name ?: key]
        def svDir = value.sv_dir ?: "${projectDir}/sv_dir/${meta.sample_name}"
        [meta, svDir]
    }
}

def xgbSamples = parseXgbSamples(params.samples, projectDir)

workflow XGBOOST_SCORING_WES {

    ch_samples = Channel.fromList(xgbSamples)
        .map { meta, svDir -> tuple(meta, file(svDir)) }

    ch_xgb_model = Channel.fromPath(
        params.xgb_model ?: "${projectDir}/modules/xgboost_consensus_model/xgboost_tuned_no_freec_cnmops.json",
        checkIfExists: true
    ).collect()

    // The one place this pipeline points at the shared modules/consensus_model/
    // (repo root) -- reused unchanged for its caller-parsing/candidate-merging
    // logic, staged into the scoring task so score_samples_xgb.py can import it
    // without hardcoding this path itself.
    ch_cnv_consensus_model_py = Channel.fromPath(
        "${projectDir}/../../modules/consensus_model/cnv_consensus_model.py",
        checkIfExists: true
    ).collect()

    XGBOOST_CONSENSUS_SCORE(ch_samples, ch_xgb_model, ch_cnv_consensus_model_py)

    CLASSIFY_AND_ANNOTATE(
        XGBOOST_CONSENSUS_SCORE.out.merged_tsv,
        ch_organism_gtf_tsv
    )

    ch_for_final = XGBOOST_CONSENSUS_SCORE.out.merged_tsv
        .map { meta, tsv, bed -> tsv }
        .mix(
            CLASSIFY_AND_ANNOTATE.out.annotated_tsv
                .map { meta, tsv -> tsv }
        )
        .collect()

    // No BAM-derived coverage step in this pipeline (input is pre-computed
    // caller calls, not BAMs) -- FINAL_FORMATTING_TABLES tolerates an empty
    // coverage list (its `coverage_out` emit is `optional: true`, and its
    // script only opportunistically globs *.region_coverage.tsv files).
    FINAL_FORMATTING_TABLES(
        ch_for_final,
        ch_organism_dna_panel,
        ch_organism_gtf,
        Channel.empty().collect()
    )
}
