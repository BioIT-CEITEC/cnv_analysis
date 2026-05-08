include { READ_DEPTH_ECOLE } from "../../modules/ECOLE/read_depth/main.nf"
include { CNV_CALL_ECOLE } from "../../modules/ECOLE/cnv_call/main.nf"

workflow ECOLE_ANALYSIS {

    take:
    ch_input_bams
    ch_regions_of_interest

    main:

    ch_bed = ch_regions_of_interest.first()

    ch_input_bams
        .map { sample ->
            sample[1..-1]
        }
        .flatten()
        .collect()
        .set { ch_all_bam_files }

    READ_DEPTH_ECOLE (
        ch_all_bam_files,
        ch_bed
    )

    ch_read_depths     = READ_DEPTH_ECOLE.out.ecole_preprocessed_samples

    CNV_CALL_ECOLE (
        ch_read_depths,
        ch_bed
    )

    emit:
    ch_ecole_varcalls = CNV_CALL_ECOLE.out.ecole_varcalls

}