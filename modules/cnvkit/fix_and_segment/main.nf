process FIX_AND_SEGMENT {

    tag "${meta.id}"

    input:
    tuple val(meta), path(tumor_target), path(tumor_antitarget)
    path cnv_reference

    output:
    tuple val(meta), path("fixed_cov.cnr"), path("segmented_cov.cns")

    script:
    """
    mkdir -p segmented
    cnvkit.py fix ${tumor_target} \
        ${tumor_antitarget} \
        ${cnv_reference} \
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