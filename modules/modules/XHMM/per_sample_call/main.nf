process PER_SAMPLE_CALL_XHMM {
    tag "${meta.sample_name}"
    publishDir "structural_varcalls/${meta.sample_name}/XHMM", mode: 'copy'
    conda "${moduleDir}/../env.yaml"

    input:
    tuple  val(meta), path(bam), path(bai)
    path cnv_calls

    output:
    tuple val(meta), path("${meta.sample_name}_xhmm.tsv"), emit: xhmm_per_sample

    script:
    """
    python3 ${projectDir}/bin/xhmm_extract_sample.py \
        --xcnv ${cnv_calls} \
        --sample-name ${bam.simpleName} \
        --output ${meta.sample_name}_xhmm.tsv
    """

    stub:
    """
    touch ${meta.sample_name}_xhmm.tsv
    """
}
