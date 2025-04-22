process DIAGRAM_AND_SCATTER_CNVKIT {

    tag "${meta.id}"
    publishDir "structural_varcalls/${meta.id}/cnvkit", mode: 'copy'

    conda "../${moduleDir}/env.yaml" 

    input:
    tuple val(meta), path(fixed_cov), path(vcf_file), path(cnvkit_cnv_calls) // mandatory: [ [meta],[cnr file],[vcf file],[cnv calls file] ]

    output:
    tuple val(meta), path("plots/")

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