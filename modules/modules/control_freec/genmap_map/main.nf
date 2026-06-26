process GENMAP_MAP_FREEC {
    conda "${moduleDir}/../env.yaml"

    input:
    path genmap_index

    output:
    path("mappability.bedgraph"), emit: mappability_bg

    script:
    """
    tmpdir=genmap_raw
    mkdir -p "\$tmpdir"
    genmap map -K ${params.freec_mappability_k} -E ${params.freec_mappability_e} \\
        -I ${genmap_index} -O "\$tmpdir" -bg \\
        -T ${task.cpus} \\
    find "\$tmpdir" -name "*.bedgraph" | sort | xargs cat \\
        > mappability.bedgraph
    rm -rf "\$tmpdir"
    """

    stub:
    """
    touch mappability.bedgraph
    """
}
