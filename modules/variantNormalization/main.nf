process VARIANT_NORMALIZATION {

    tag "${sample_id}"
    publishDir "variantNormalization/${sample_id}", mode: 'copy'
    conda "${moduleDir}/env.yaml"

    input:
    tuple val(sample_id), path(input_files)

    output:
    tuple val(sample_id), path("${sample_id}_normalized"), emit: normalized_varcalls

    script:
    """
    python ${projectDir}/bin/variant_normalization.py \
        ./ ${sample_id}_normalized
    """
}
