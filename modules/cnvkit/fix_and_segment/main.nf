process FIX_AND_SEGMENT_CNVKIT {

    tag "${meta.donor}"

    conda "${moduleDir}/../env.yaml"

    input:
    tuple val(meta), path(tumor_target), path(tumor_antitarget), path(normal_target), path(normal_antitarget) // mandatory [ [meta],[target],[antitarget] ] target and antitarget files produced in the coverage process
    path cnvkit_reference // mandatory path to the reference file "normal_reference.cnn" produced at reference process

    output:
    tuple val(meta), path("segmented/fixed_cov.cnr"), path("segmented/segmented_cov.cns"), emit: cnvkit_segments

    script:
    """
    mkdir -p segmented
    cnvkit.py fix ${tumor_target} \
        ${tumor_antitarget} \
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