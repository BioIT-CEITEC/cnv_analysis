process VARIANT_NORMALIZATION_JABCONTOOL {

    publishDir "structural_varcalls/all_samples/", mode: 'copy'
    conda "${moduleDir}/env.yaml"

    input:
    path(input_dir)

    output:
    path("jabCoNtool/normalized_varcalls/*.tsv"), emit: normalized_varcalls

    script:
    """
    python ${projectDir}/bin/variant_normalization.py \
        ./ jabCoNtool/normalized_varcalls
    """

    stub:
    """
    mkdir -p jabCoNtool/normalized_varcalls
    touch jabCoNtool/normalized_varcalls/placeholder.tsv
    """
}