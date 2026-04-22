include { READ_DEPTH_ECOLE } from "../../modules/ECOLE/read_depth/main.nf"
include { CNV_CALL_ECOLE } from "../../modules/ECOLE/cnv_call/main.nf"

workflow ECOLE_ANALYSIS {

    take:
    ch_input_bams
    ch_regions_of_interest

    main:


    // Transform input channel
    ch_input_bams
        .map { sample ->
            sample[1..-1] // Remove the first element if necessary
        }
        .flatten() // Flatten the nested structure
        .collect() // Collect all files into a single list
        .set { ch_all_bam_files }

    // Debug transformed channel

    READ_DEPTH_ECOLE (
        ch_all_bam_files,
        ch_regions_of_interest
    )

    ch_read_depths     = READ_DEPTH_ECOLE.out.ecole_preprocessed_samples

    CNV_CALL_ECOLE (
        ch_read_depths,
        ch_regions_of_interest
    )

}