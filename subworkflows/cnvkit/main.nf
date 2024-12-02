include Utils

import { PREPARE_REGIONS_CNVKIT as PREPARE_REGIONS } from "../../modules/cnvkit/prepare_regions/main" 
import { GET_COVERAGE_CNVKIT as GET_COVERAGE } from "../../modules/cnvkit/coverage/main"
import { FIX_AND_SEGMENT_CNVKIT as FIX_AND_SEGMENT } from "../../modules/cnvkit/fix_and_segment/main"
import { CNV_CALL_CNVKIT as CNV_CALL } from "../../modules/cnvkit/call/main"
import { DIAGRAM_AND_SCATTER_CNVKIT as DIAGRAM_AND_SCATTER } from "../../modules/cnvkit/plotting/main"

workflow CNKIT {

    take:
    // Inputs from config

   ch_inputs // [mandatory] [ meta, [tumor_bam], [tumor_bam_bai], [normal_bam], [normal_bam_bai] ]




 
}