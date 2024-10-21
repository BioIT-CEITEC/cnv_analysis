nextflow.enable.dsl = 2


// Check that all the references and input files exist

if (file(params.input).exists()) { ch_input_data = file(params.input) } else { exit 1, 'Input samplesheet not found in the directory!'}
if (file(params.organism_fasta).exists()) { ch_organism_fasta = file(params.organism_fasta) } else { exit 1, 'The reference genome does not exist'}
if (file("${params.organism_fasta}.fai").exists()) { ch_organism_fai = file(${params.organism_fasta}.fai) } else { exit 1, 'The reference genome index does not exist'}
if (file(params.organism_dict).exists()) { ch_organism_dict = file(params.organism_dict) } else { exit 1, 'The reference genome dictionary does not exist'}

include { INPUT_CHECK } from './subworkflows/input_check.nf'
include { PREPROCESSING } from './modules/preprocess.nf'
include { COHORT_PREPROCESS } from './modules/cohort_data.nf'
include { PURPLE_ANALYSIS } from './workflows/purple.nf'
include { CONTROL_FREEC } from './modules/control_freec.nf'
include { JABCONTOOL_ANALYSIS } from './workflows/jabCoNtool.nf'


workflow CNV_ANALYSIS {

    INPUT_CHECK (
        ch_input_data
    )

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