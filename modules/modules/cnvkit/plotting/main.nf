process DIAGRAM_AND_SCATTER_CNVKIT {

    tag "${meta.sample_name}"
    publishDir "structural_varcalls/${meta.sample_name}", mode: 'copy'

    conda "${moduleDir}/../env.yaml"

    input:
    path(cnvkit_cnv_calls)
    tuple path(fixed_cov), path(segmented_cov)

    output:
    tuple val(meta), path("cnvkit/*")

    script:

    """
    mkdir -p cnvkit
    cnvkit.py diagram ${fixed_cov} -s ${cnvkit_cnv_calls} -o cnvkit/cnvkit_diagram.pdf
    cnvkit.py scatter ${fixed_cov} -s ${cnvkit_cnv_calls} -v ${vcf_file} -o cnvkit/cnvkit_scatter.png
    """

    stub:
    """
    mkdir -p cnvkit
    touch cnvkit/cnvkit_diagram.pdf
    touch cnvkit/cnvkit_scatter.png
    """

}