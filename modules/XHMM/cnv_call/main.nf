process CNV_CALL_XHMM {
    publishDir "structural_varcalls/all_samples", mode: 'copy'
    conda "${moduleDir}/../env.yaml"

    input:
    path bams
    path target_bed
    val sample_names

    output:
    path("xhmm_calls_output/*.xcnv"),    emit: xhmm_cnvcalls
    path("xhmm_calls_output/xhmm_work"), emit: xhmm_intermediates

    script:
    def bam_files_str = bams.findAll { it.name.endsWith('.bam') }.join(' ')
    def names_str     = (sample_names instanceof List ? sample_names : [sample_names]).join(',')
    def target_bed_path = target_bed instanceof Collection ? target_bed[0] : target_bed
    def raw_params_values = params.xhmm_params_values
    def params_values = raw_params_values instanceof List ? raw_params_values.collect { it.toString() } : raw_params_values.toString()
        .replaceAll(/[\[\]]/, '')
        .split(/\s*,\s*|\s+/)
        .findAll { it }
    """
    mkdir -p xhmm_calls_output

    python3 ${projectDir}/bin/xhmm_cohort_call.py \
        --bams ${bam_files_str} \
        --sample-names ${names_str} \
        --bed ${target_bed} \
        --output xhmm_calls_output/calls.xcnv \
        --threads 1 \
        --min-target-size ${params.xhmm_min_target_size} \
        --max-target-size ${params.xhmm_max_target_size} \
        --min-mean-target-rd ${params.xhmm_min_mean_target_rd} \
        --max-mean-target-rd ${params.xhmm_max_mean_target_rd} \
        --min-mean-sample-rd ${params.xhmm_min_mean_sample_rd} \
        --max-mean-sample-rd ${params.xhmm_max_mean_sample_rd} \
        --max-sd-sample-rd ${params.xhmm_max_sd_sample_rd} \
        --pve-mean-factor ${params.xhmm_pve_mean_factor} \
        --max-sd-target-rd ${params.xhmm_max_sd_target_rd} \
        --discover-some-qual-thresh ${params.xhmm_discover_some_qual_threshold} \
        --params-values ${params_values.join(' ')}
    """

    stub:
    """
    mkdir -p xhmm_calls_output/xhmm_work
    touch xhmm_calls_output/calls.xcnv
    """
}
