process CNVKIT_CALL {

    tag "${meta.donor}"
    publishDir "structural_varcalls/${meta.donor}/cnvkit", mode: 'copy'

    conda "${moduleDir}/../env.yaml" 

    input:
    tuple val(meta), path(fixed_cov), path(segmented_cov)

    output:
    path("call/CNV_calls.cns"), emit: cnvkit_calls

    script:

    """
    mkdir -p call
    cnvkit.py call -y -m clonal ${segmented_cov} -o call/CNV_calls.cns --purity 0.5
    """

    stub:
    """
    mkdir -p call
    touch call/CNV_calls.cns
    """
}