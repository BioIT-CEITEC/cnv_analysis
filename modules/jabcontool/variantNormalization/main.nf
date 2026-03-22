process VARIANT_NORMALIZATION_JABCONTOOL {

    publishDir "variantNormalization/jabcontool", mode: 'copy'
    conda "${moduleDir}/env.yaml"

    input:
    path(input_dir)

    output:
    path("normalized/*"), emit: normalized_varcalls

    script:
    """
    python ${projectDir}/bin/variant_normalization.py \
        ./ normalized
    """
}