process FIX_AND_SEGMENT_CNVKIT {

    tag "${meta.sample_name}"

    conda "${moduleDir}/../env.yaml"

    input:
    tuple val(meta), path(target), path(antitarget)
    path cnvkit_reference // mandatory path to the reference file "normal_reference.cnn" produced at reference process

    output:
    tuple val(meta), path("segmented/fixed_cov.cnr"), path("segmented/segmented_cov.cns"), emit: cnvkit_segments

    script:
    """
    mkdir -p segmented
    cnvkit.py fix ${target} \
        ${antitarget} \
        ${cnvkit_reference} \
        -o segmented/fixed_cov.cnr

    cnvkit.py segment segmented/fixed_cov.cnr \
        -o segmented/segmented_cov.cns
    """

    stub:
    """
    mkdir -p segmented
    touch segmented/fixed_cov.cnr
    touch segmented/segmented_cov.cns
    """
}