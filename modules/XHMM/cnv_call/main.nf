process CNV_CALL_XHMM {
    conda "${moduleDir}/../env.yaml"

    input:
    path bams
    path target_bed

    output:
    path("xhmm_calls_output/*.xcnv"), emit: xhmm_cnvcalls

    script:
    def bam_files_str = bams.findAll { it.name.endsWith('.bam') }.join(' ')
    def names_str     = (sample_names instanceof List ? sample_names : [sample_names]).join(',')
    """
    mkdir -p xhmm_calls_output

    python3 ${projectDir}/bin/xhmm_cohort_call.py \
        --bams ${bam_files_str} \
        --sample-names ${names_str} \
        --bed ${target_bed} \
        --output xhmm_calls_output/calls.xcnv \
        --threads ${task.cpus}
        #--min-target-size ${params.xhmm_min_target_size}
        #--max-target-size ${params.xhmm_max_target_size}
        #--min-mean-target-rd ${params.xhmm_min_mean_target_rd}
        #--max-mean-target-rd ${params.xhmm_max_mean_target_rd}
        #--min-mean-sample-rd ${params.xhmm_min_mean_sample_rd}
        #--max-mean-sample-rd ${params.xhmm_max_mean_sample_rd}
        #--max-sd-sample-rd ${params.xhmm_max_sd_sample_rd}
        #--pve-mean-factor ${params.xhmm_pve_mean_factor}
        #--max-sd-target-rd ${params.xhmm_max_sd_target_rd}
        #--discover-some-qual-threshold ${params.xhmm_discover_some_qual_threshold}
        #--params-values ${params.xhmm_params_values}
    """

    stub:
    """
    mkdir -p xhmm_calls_output
    touch xhmm_calls_output/calls.xcnv
    """
}
