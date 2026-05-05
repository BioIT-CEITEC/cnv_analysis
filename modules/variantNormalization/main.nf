process VARIANT_NORMALIZATION {

    tag "${meta.sample_name}"
    publishDir "structural_varcalls/${meta.sample_name}", mode: 'copy'
    conda "${moduleDir}/env.yaml"
    cache false

    input:
    tuple val(meta), path(input_dir)

    output:
    tuple val(meta), path("normalized_varcalls/${meta.sample_name}_*_normalized.tsv"), emit: normalized_varcalls

    script:
    """
    mkdir -p normalized_varcalls
    python ${projectDir}/bin/variant_normalization.py \
        ./ normalized_varcalls
    """

    stub:
    """
    mkdir -p normalized_varcalls
    touch normalized_varcalls/${meta.sample_name}_normalized.tsv
    """
}