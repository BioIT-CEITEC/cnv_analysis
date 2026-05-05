process CNVKIT_CALL {

    tag "${meta.sample_name}"
    publishDir "structural_varcalls/${meta.sample_name}", mode: 'copy'

    conda "${moduleDir}/../env.yaml" 

    input:
    tuple val(meta), path(fixed_cov), path(segmented_cov)

    output:
    tuple val(meta), path("cnvkit/*.cns"), path("cnvkit/${meta.sample_name}_preannot_calls.bed"), path("cnvkit/${meta.sample_name}_CNV_calls.bed"), path("cnvkit/*.vcf"), path("cnvkit/${meta.sample_name}_cnvkit.tsv"), emit: cnvkit_calls

    script:

    """
    mkdir -p cnvkit
    cnvkit.py call ${segmented_cov} -m threshold --purity 0.5 -o cnvkit/${meta.sample_name}_CNV_calls.cns
    cnvkit.py export bed cnvkit/${meta.sample_name}_CNV_calls.cns --show all -o cnvkit/${meta.sample_name}_CNV_calls.bed
    cnvkit.py export vcf cnvkit/${meta.sample_name}_CNV_calls.cns -i ${meta.sample_name} -o cnvkit/${meta.sample_name}_CNV_calls.vcf
    vcf2tsv.py cnvkit/${meta.sample_name}_CNV_calls.vcf cnvkit/${meta.sample_name}_cnvkit.tsv cnvkit/${meta.sample_name}_preannot_calls.bed
    """

    stub:
    """
    mkdir -p cnvkit
    touch cnvkit/${meta.sample_name}_CNV_calls.cns
    """
}