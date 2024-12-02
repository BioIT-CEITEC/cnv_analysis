#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here's very important where major changes need to me made if at some
point we decide to try to submit it to nf-core. It has its specific
built-in functions that are not being used here as the way of loading
data, references and config differs.
*/

include { TARGETED } from './workflows/targeted'
include { WGS     } from './workflows/wgs'

/*
-----------------------------------------------------------------------
    WORKFLOW DEFINITION
-----------------------------------------------------------------------
*/

workflow CNV_ANALYSIS {

    if (params.lib_ROI === "wgs") {
        WGS()
    } else {
        TARGETED()
    }
}

/*
-----------------------------------------------------------------------
    PIPELINE EXECUTION
-----------------------------------------------------------------------
*/

workflow {
    CNV_ANALYSIS()
}
