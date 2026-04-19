process DIFF_REFERENCE {
    tag "${gene_region}_${pseudogene_region}"
    publishDir "pseudogene/test", mode: 'copy'
    conda "${moduleDir}/../env.yaml"

    input:
    tuple val(gene_region), val(pseudogene_region), path(gene_ref), path(pseudogene_ref)

    output:
    tuple val(gene_region), val(pseudogene_region), path("${gene_ref.baseName}_${pseudogene_ref.baseName}_diff.tsv"), emit: diff_output

    script:
    def use_global = params.use_global == null ? true : params.use_global
    """
    if [ "${use_global}" == "true" ]; then
        python ${projectDir}/bin/diff_reference.py --gene ${gene_ref} --pseudogene ${pseudogene_ref} --output "${gene_ref.baseName}_${pseudogene_ref.baseName}_diff.tsv" --use_global
    else
        python ${projectDir}/bin/diff_reference.py --gene ${gene_ref} --pseudogene ${pseudogene_ref} --output "${gene_ref.baseName}_${pseudogene_ref.baseName}_diff.tsv"
    fi
    """
}