#!/usr/bin/env nextflow

process VARIANT_NORMALIZATION {

    publishDir "variantNormalization", mode: 'copy'
    conda "${moduleDir}/env.yaml"

    input:
    tuple val(meta), path(input_dir)

    output:
    tuple val(meta), path("normalized"), emit: normalized

    script:
    """
    python ${projectDir}/bin/variant_normalization.py \
        ${input_dir} normalized
    """
}
