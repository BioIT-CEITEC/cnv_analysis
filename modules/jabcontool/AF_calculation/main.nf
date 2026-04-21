process SNP_AF_CALC {

    tag "${meta.sample_name}"
    conda "${moduleDir}/env.yaml"

    input:
    tuple val(meta), path(bam), path(bai)
    path reference_fasta
    path reference_index
    path reference_snps

    output:
    tuple val(meta), path("${meta.sample_name}.snpAF.tsv"), emit: snpAF

    script:

    """
    alleleCounter -r ${reference_fasta} -l ${reference_snps} -b ${bam} -o ${meta.sample_name}.snpAF.tsv
    """
    stub:
    """
    touch ${meta.sample_name}.snpAF.tsv
    """
}