process CNVKIT_CALL {

    tag "${meta.id}"
    publishDir "structural_varcalls/${meta.id}/cnvkit", mode: 'copy'

    conda "../${moduleDir}/env.yaml" 

    input:
    tuple val(meta), path(vcf_file), path(segmented_cov)

    output:
    tuple val(meta), path("call/")

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