process DiffReference {
    tag "${gene_region}_${pseudogene_region}"
    conda './env.yml'
    input:
    tuple val(gene_region), val(pseudogene_region), path(gene_ref), path(pseudogene_ref)

    output:
    tuple val(gene_region), val(pseudogene_region), path("${gene_ref.baseName}_${pseudogene_ref.baseName}_diff.tsv")

    script:
    """
    if [ "${params.use_global}" == "true" ]; then
        diff_reference.py --gene ${gene_ref} --pseudogene ${pseudogene_ref} --output "${gene_ref.baseName}_${pseudogene_ref.baseName}_diff.tsv" --use_global
    else
        diff_reference.py --gene ${gene_ref} --pseudogene ${pseudogene_ref} --output "${gene_ref.baseName}_${pseudogene_ref.baseName}_diff.tsv"
    fi
    """
}