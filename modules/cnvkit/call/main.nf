process CNVKIT_CALL {

    tag "${meta.sample_name}"
    publishDir "structural_varcalls/cnvkit", mode: 'copy'

    conda "${moduleDir}/../env.yaml" 

    input:
    tuple val(meta), path(fixed_cov), path(segmented_cov)

    output:
    tuple val(meta), path("call/*.cns"), path("call/${meta.sample_name}_preannot_calls.bed"), path("call/${meta.sample_name}_CNV_calls.bed"), path("call/*.vcf"), path("call/*.tsv"), emit: cnvkit_calls

    script:

    """
    mkdir -p call
    cnvkit.py call ${segmented_cov} -m threshold --purity 0.5 -o call/${meta.sample_name}_CNV_calls.cns
    cnvkit.py export bed call/${meta.sample_name}_CNV_calls.cns --show all -o call/${meta.sample_name}_CNV_calls.bed
    cnvkit.py export vcf call/${meta.sample_name}_CNV_calls.cns -i ${meta.sample_name} -o call/${meta.sample_name}_CNV_calls.vcf
    vcf2tsv.py call/${meta.sample_name}_CNV_calls.vcf call/${meta.sample_name}_preannot_calls.tsv call/${meta.sample_name}_preannot_calls.bed
    """

    stub:
    """
    mkdir -p call
    touch call/CNV_calls.cns
    """
}