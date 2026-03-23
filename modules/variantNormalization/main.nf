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

        manifest=${meta.sample_name}_normalized/produced_outputs.tsv
        printf "tool\tpath\n" > "$manifest"

        [[ -f ${meta.sample_name}_normalized/${meta.sample_name}_cnMOPS_normalized.tsv ]] && \
            printf "cnMOPS\t%s\n" "${meta.sample_name}_normalized/${meta.sample_name}_cnMOPS_normalized.tsv" >> "$manifest"
        [[ -f ${meta.sample_name}_normalized/${meta.sample_name}_cnvkit_normalized.tsv ]] && \
            printf "cnvkit\t%s\n" "${meta.sample_name}_normalized/${meta.sample_name}_cnvkit_normalized.tsv" >> "$manifest"
        [[ -f ${meta.sample_name}_normalized/${meta.sample_name}_ExomeDepth_normalized.tsv ]] && \
            printf "ExomeDepth\t%s\n" "${meta.sample_name}_normalized/${meta.sample_name}_ExomeDepth_normalized.tsv" >> "$manifest"
        [[ -f ${meta.sample_name}_normalized/${meta.sample_name}_panelcnMOPS_normalized.tsv ]] && \
            printf "panelcnMOPS\t%s\n" "${meta.sample_name}_normalized/${meta.sample_name}_panelcnMOPS_normalized.tsv" >> "$manifest"
        [[ -f ${meta.sample_name}_normalized/${meta.sample_name}_gatk_normalized.tsv ]] && \
            printf "gatk\t%s\n" "${meta.sample_name}_normalized/${meta.sample_name}_gatk_normalized.tsv" >> "$manifest"
    """
}