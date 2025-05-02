nextflow.enable.dsl = 2


// Validate inputs and potentially references

include { CONTROL_FREEC } from "../modules/control_freec/main.nf"
include { BED_PREPARATION } from "../modules/preprocessing/bed_preparation/main.nf"
include { COHORT_PREPARATION } from "../modules/preprocessing/cohort_preparation/main.nf"
include { CNVKIT_ANALYSIS } from "../subworkflows/cnvkit/main.nf"
include { JABCONTOOL_ANALYSIS } from "../subworkflows/jabcontool/main.nf"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REFERENCE CONFIG LOADING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Load the reference data from the paths in config file
params = Utils.load_organism(params)

ch_organism_fasta = params.organism_fasta ? Channel.fromPath(params.organism_fasta).collect() : Channel.empty()
ch_organism_fasta_fai = params.organism_fasta ? Channel.fromPath(params.organism_fasta + '.fai').collect() : Channel.empty()
organism_dna_panel = params.organism_dna_panel ? Channel.fromPath(params.organism_dna_panel).collect() : Channel.empty()
organism_cytoband = params.organism_cytoband ? Channel.fromPath(params.organism_cytoband).collect() : Channel.empty()

inputs = Utils.parseInputVC(params.input, params.normal_tumor, log)

workflow TARGETED {

    // Create the channel from the parseInputVC function
    // channel: [ meta, []]
    ch_inputs = Channel.fromList(inputs)

    BED_PREPARATION (
        ch_organism_fasta,
        ch_organism_fasta_fai
    )

    ch_binned_genome = BED_PREPARATION.out.binned_genome
    ch_gc_profile = BED_PREPARATION.out.gc_profile

/*
    if (!cohort_data_channel.empty) {
        COHORT_PREPROCESS (
            cohort_data_channel
        )
    }

    if (params.lib_ROI == "wgs") {

        PURPLE_ANALYSIS (
            ch_input_bams
        )

        ch_purple_results = PURPLE_ANALYSIS.out.ch_purple_outputs

        CONTROL_FREEC (
            ch_input_bams,
            ch_wgs_preprocess
        )
        ch_control_freec_results = CONTROL_FREEC.out.var_call

    } else {

        CNVKIT_ANALYSIS (
            ch_input_bams
        )

    }
*/
}