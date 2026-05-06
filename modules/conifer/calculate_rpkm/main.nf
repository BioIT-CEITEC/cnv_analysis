process CALCULATE_RPKM_CONIFER {
    tag "${meta.sample_name}"
    conda "${moduleDir}/../env.yaml"

    input:
    tuple val(meta), path(bam), path(bai)
    path probes

    output:
    tuple val(meta), path("${meta.sample_name}.rpkm.txt"), emit: rpkm_conifer

    script:
    def conifer = "${moduleDir}/../conifer_v0.2.2/conifer.py"
    """
    python ${conifer} \
        rpkm \
        --probes ${probes} \
        --input ${bam} \
        --output ${meta.sample_name}.rpkm.txt
    """

    stub:
    """
    touch ${meta.sample_name}.rpkm.txt
    """
}
