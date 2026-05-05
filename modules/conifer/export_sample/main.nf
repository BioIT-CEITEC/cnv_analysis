process EXPORT_SAMPLE_CONIFER {
    tag "${meta.sample_name}"
    publishDir "structural_varcalls/${meta.sample_name}/conifer", mode: 'copy'
    conda "${moduleDir}/../env.yaml"

    input:
    tuple val(meta), path(cohort_hdf5)

    output:
    tuple val(meta), path("${meta.sample_name}_conifer.bed"), emit: conifer_export

    script:
    def conifer = "${moduleDir}/../conifer_v0.2.2/conifer.py"
    def out     = "${meta.sample_name}_conifer.bed"
    """
    exported=0
    for candidate in "${meta.sample_name}" "${meta.sample_name}.rpkm"; do
        python ${conifer} export \
            --input  ${cohort_hdf5} \
            --sample "\$candidate" \
            --output ${out} && exported=1 && break
    done
    [ "\$exported" -eq 1 ] || touch ${out}
    """

    stub:
    """
    touch ${meta.sample_name}_conifer.bed
    """
}
