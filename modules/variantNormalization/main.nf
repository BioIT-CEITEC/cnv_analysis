process VARIANT_NORMALIZATION {

    tag "${meta.sample_name}"
    publishDir "variantNormalization/", mode: 'copy'
    conda "${moduleDir}/env.yaml"

    input:
    tuple val(meta), path(input_dir)

    output:
    tuple val(meta), path("${meta.sample_name}_normalized/${meta.sample_name}_*_normalized.tsv"), emit: normalized_varcalls

    script:
    """
    mkdir -p ${meta.sample_name}_normalized
    python ${projectDir}/bin/variant_normalization.py \
        ./ ${meta.sample_name}_normalized
    """
}