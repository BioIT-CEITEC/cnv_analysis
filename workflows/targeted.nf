nextflow.enable.dsl = 2

// Check if the specified input samplesheet exists
if (file(params.input).exists()) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not found in the directory!'}

// Get input information from samplesheet
inputs = Utils.processInput(params.input)

// Validate inputs and potentially references

include { AMBER_ANALYSIS } from './subworkflows/amber_analysis.nf'
include { COBALT_ANALYSIS } from './subworkflows/cobalt_analysis.nf'
include { GRIDSS_ANALYSIS } from './subworkflows/gridss_analysis.nf'
include { GRIPSS_ANALYSIS } from './subworkflows/gripss_analysis.nf'
include { PURPLE_CALL } from './subworkflows/purple_call.nf'

workflow WGS {

    // Create the channel from the samplesheet CSV
    // channel: [ meta, []]
    ch_inputs = Channel.fromList(inputs)

    AMBER_ANALYSIS(
        ch_inputs
    )




}




workflow CNV_ANALYSIS {

    AMBER_ANALYSIS {
        ch_input_bams
    }

    ch_input_bams = INPUT_CHECK.out.samples

    PREPROCESSING {
        ch_organism_fasta
        ch_organism_fai
    }

    def ch_cohort_input = params.use_cohort_data ? file(params.cohort_data) : Channel.empty()

    if (params.use_cohort_data) {
        COHORT_PREPROCESS (
            ch_cohort_input
        )
    }

    def ch_wgs_preprocess = params.lib_ROI == "wgs" ? PREPROCESSING.out.preprocessed : Channel.empty()

    if (params.lib_ROI == "wgs") {

        PURPLE_ANALYSIS {
            ch_input_bams
        }

        ch_purple_results = PURPLE_ANALYSIS.out.ch_purple_outputs

        CONTROL_FREEC {
            ch_input_bams
            ch_wgs_preprocess
        }
        ch_control_freec_results = CONTROL_FREEC.out.var_call

    } else {

        CNVKIT_ANALYSIS (
            ch_input_bams
        )



    }

    def ch_cohort_data = params.use_cohort_data ? COHORT_PREPROCESS.out.cohort_data : Channel.empty()

    JABCONTOOL_ANALYSIS {
        ch_input_bams
        ch_wgs_preprocess
        ch_cohort_data
    }



}