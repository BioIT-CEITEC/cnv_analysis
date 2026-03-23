process VARIANT_NORMALIZATION {

    tag "${meta.sample_name}"
    publishDir "variantNormalization/", mode: 'copy'
    conda "${moduleDir}/env.yaml"

    input:
    tuple val(meta), path(input_dir)

    output:
    tuple val(meta), 
    path("${meta.sample_name}_normalized/${meta.sample_name}_cnMOPS_normalized.tsv", optional: true),
    path("${meta.sample_name}_normalized/${meta.sample_name}_cnvkit_normalized.tsv", optional: true),
    path("${meta.sample_name}_normalized/${meta.sample_name}_ExomeDepth_normalized.tsv", optional: true),
    path("${meta.sample_name}_normalized/${meta.sample_name}_panelcnMOPS_normalized.tsv", optional: true),
    path("${meta.sample_name}_normalized/${meta.sample_name}_gatk_normalized.tsv", optional: true), emit: normalized_varcalls

    script:
    """
    python ${projectDir}/bin/variant_normalization.py \
        ./ ${meta.sample_name}_normalized
    """
}