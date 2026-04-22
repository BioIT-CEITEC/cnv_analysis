process EXTRACT_SAMPLE_CONIFER {
    tag "${meta.sample_name}"
    publishDir "${params.publish_dir}/conifer", mode: 'copy', overwrite: true
    conda "${moduleDir}/../env.yaml"

    input:
    tuple val(meta), path(bam), path(bai)
    path  cohort_calls

    output:
    tuple val(meta), path("${meta.sample_name}_conifer.tsv"), emit: conifer_per_sample

    script:
    """
    python3 ${projectDir}/bin/conifer_extract_sample.py \
        --calls       ${cohort_calls} \
        --sample-name ${meta.sample_name} \
        --output      ${meta.sample_name}_conifer.tsv
    """

    stub:
    """
    touch ${meta.sample_name}_conifer.tsv
    """
}
