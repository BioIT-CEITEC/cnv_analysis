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
Utils.set_data_tags(params)
Utils.loadSample(params)

ch_organism_fasta = params.organism_fasta ? Channel.fromPath(params.organism_fasta).collect() : Channel.empty()
ch_organism_fasta_fai = params.organism_fasta ? Channel.fromPath(params.organism_fasta + '.fai').collect() : Channel.empty()
ch_organism_dict = params.organism_dict ? Channel.fromPath(params.organism_dict).collect() : Channel.empty()
ch_organism_dna_panel = params.organism_dna_panel ? Channel.fromPath(params.organism_dna_panel).collect() : Channel.empty()
ch_organism_cytoband = params.organism_cytoband ? Channel.fromPath(params.organism_cytoband).collect() : Channel.empty()
ch_organism_snps = params.organism_snps_panel ? Channel.fromPath(params.organism_snps_panel).collect() : Channel.empty()



inputs = Utils.parseInputVC(params.new_samples, params.normal_tumor, projectDir, log)


workflow TARGETED {

    ch_cohort_data = Channel.empty()

    // Create the channel from the parseInputVC function
    // channel: [ meta, []]
    ch_inputs = Channel.fromList(inputs)

    BED_PREPARATION (
        ch_organism_fasta,
        ch_organism_dict
    )

    ch_binned_genome = BED_PREPARATION.out.binned_genome
    ch_gc_profile = BED_PREPARATION.out.gc_profile

    CNVKIT_ANALYSIS (

    ch_inputs,
    ch_organism_dna_panel,
    ch_organism_fasta,
    ch_organism_fasta_fai

    )

    JABCONTOOL_ANALYSIS (
    ch_inputs,
    ch_organism_fasta,
    ch_organism_fasta_fai,
    ch_organism_snps,
    ch_organism_dna_panel,
    ch_organism_cytoband,
    ch_gc_profile,
    ch_cohort_data

    )

    CONTROL_FREEC (
    ch_inputs,
    ch_organism_fasta,
    ch_organism_fasta_fai,
    ch_organism_snps,
    ch_binned_genome,
    ch_gc_profile

    )

}