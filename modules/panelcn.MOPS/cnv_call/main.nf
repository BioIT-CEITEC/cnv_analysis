process CNV_CALL_PANELCNMOPS {

    tag "${meta.sample_name}"
    publishDir "structural_varcalls/${meta.sample_name}", mode: 'copy'
    conda "${moduleDir}/env.yaml"

    input:
    tuple val(meta), path(bam), path(bai)
    path cohort_data

    output:
    tuple val(meta), path("panelcnMOPS/${meta.sample_name}_panelcnMOPS.tsv"), emit: panelcnMOPS_cnvcalls

    script:
    """
    mkdir -p panelcnMOPS
    Rscript ${projectDir}/bin/panelcnMOPS_wrapper.R \
    ${bam} \
    ${cohort_data} \
    ${meta.sample_name} \
    panelcnMOPS/${meta.sample_name}_panelcnMOPS.tsv
    """

    stub:
    """
    mkdir -p panelcnMOPS
    touch panelcnMOPS/${meta.sample_name}_panelcnMOPS.tsv
    """
}