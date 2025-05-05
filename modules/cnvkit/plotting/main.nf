process DIAGRAM_AND_SCATTER_CNVKIT {

    tag "${meta.donor}"
    publishDir "structural_varcalls/${meta.donor}/cnvkit", mode: 'copy'

    conda "${moduleDir}/../env.yaml"

    input:
    path(cnvkit_cnv_calls)
    tuple path(fixed_cov), path(segmented_cov)

    output:
    tuple val(meta), path("plots/*")

    script:

    """
    mkdir -p plots
    cnvkit.py diagram ${fixed_cov} -s ${cnvkit_cnv_calls} -o plots/cnvkit_diagram.pdf
    cnvkit.py scatter ${fixed_cov} -s ${cnvkit_cnv_calls} -v ${vcf_file} -o plots/cnvkit_scatter.png
    """

    stub:
    """
    mkdir -p plots
    touch plots/cnvkit_diagram.pdf
    touch plots/cnvkit_scatter.png
    """

}