process VARIANT_NORMALIZATION {

    tag "${meta.id}"
    publishDir "variantNormalization/${meta.id}", mode: 'copy'
    conda "${moduleDir}/env.yaml"

    input:
    tuple val(meta), path(input_files)

    output:
    tuple val(meta), path("${meta.id}_normalized"), emit: normalized_varcalls

    script:
    """
    python ${projectDir}/bin/variant_normalization.py \
        ./ ${meta.id}_normalized
    """
}
