process VARIANT_NORMALIZATION {

    tag "${meta.sample_name}"
    publishDir "variantNormalization/${meta.sample_name}_normalized", mode: 'copy'
    conda "${moduleDir}/env.yaml"

    input:
    tuple val(meta), path(input_dir)

    output:
    tuple val(meta), path("${meta.sample_name}_normalized/*"), emit: normalized_varcalls

    script:
    """
    python ${projectDir}/bin/variant_normalization.py \
        ./ ${meta.sample_name}_normalized
    """
}