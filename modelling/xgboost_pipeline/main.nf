#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
Standalone pipeline: score pre-computed CNV caller output files with the
tuned XGBoost consensus model (no freec/cnMOPS), then annotate and produce
the same cohort-level tables/HTML report the main CNV_ANALYSIS pipeline
does. Does not run the callers itself -- see workflows/XGBOOST_SCORING_WES.nf
and nextflow.config for the expected input layout.
*/

include { XGBOOST_SCORING_WES } from './workflows/XGBOOST_SCORING_WES.nf'

workflow XGBOOST_CNV_SCORING {
    XGBOOST_SCORING_WES()
}

workflow {
    XGBOOST_CNV_SCORING()
}
