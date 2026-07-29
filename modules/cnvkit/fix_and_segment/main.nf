process FIX_AND_SEGMENT_CNVKIT {

    tag "${meta.sample_name}"

    conda "${moduleDir}/../env.yaml"

    input:
    tuple val(meta), path(target), path(antitarget)
    path cnvkit_reference

    output:
    tuple val(meta), path("segmented/fixed_cov.cnr"), path("segmented/segmented_cov.cns"), emit: cnvkit_segments

    script:
    """
    mkdir -p segmented
    cnvkit.py fix ${target} \
        ${antitarget} \
        ${cnvkit_reference} \
        -o segmented/fixed_cov.cnr

    # Remove bins with NaN or -Inf log2 values — CBS (R/DNAcopy) crashes on them
    awk 'NR==1 || (\$5 != "NaN" && \$5 != "-Inf" && \$5 != "Inf")' segmented/fixed_cov.cnr \
        > segmented/fixed_cov_filtered.cnr

    cnvkit.py segment segmented/fixed_cov_filtered.cnr \
        --method cbs \
        --drop-low-coverage \
        -o segmented/segmented_cov_raw.cns

    awk 'NR==1 || \$7 >= 3' segmented/segmented_cov_raw.cns > segmented/segmented_cov.cns
    """

    stub:
    """
    mkdir -p segmented
    touch segmented/fixed_cov.cnr
    touch segmented/segmented_cov.cns
    """
}