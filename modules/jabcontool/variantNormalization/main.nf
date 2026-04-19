process VARIANT_NORMALIZATION_JABCONTOOL {

    publishDir "variantNormalization/", mode: 'copy'
    conda "${moduleDir}/env.yaml"

    input:
    path(input_dir)

    output:
    path("jabcontool/*.tsv"), emit: normalized_varcalls

    script:
    """
    python ${projectDir}/bin/variant_normalization.py \
        ./ jabcontool
    """
}