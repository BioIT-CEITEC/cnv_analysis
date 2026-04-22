process CALCULATE_RPKM_CONIFER {
    tag "${meta.sample_name}"
    conda "${moduleDir}/../env.yaml"

    input:
    tuple val(meta), path(bam), path(bai)
    path probes

    output:
    tuple val(meta), path("${meta.sample_name}_rpkm.txt"), emit: rpkm_conifer

    script:
    """
    python3 ${moduleDir}/conifer_v0.2.2/conifer.py \
        rpkm \
        --probes ${probes} \
        --input ${bam} \
        --output ${meta.sample_name}_rpkm.txt
    """

    stub:
    """
    touch ${meta.sample_name}_rpkm.txt
    """
}
