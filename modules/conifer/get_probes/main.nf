process GET_PROBES_CONIFER {
    conda "${moduleDir}/../env.yaml"

    input:
    path bed
    path gtf

    output:
    path "conifer_probes.tsv", emit: conifer_probes

    script:
    """
    python3 ${projectDir}/bin/conifer_prepare_probes.py \
        --bed    ${bed} \
        --gtf    ${gtf} \
        --output conifer_probes.tsv
    """

    stub:
    """
    touch conifer_probes.tsv
    """
}
