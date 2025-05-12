process SNP_AF_CALC {

    tag "${meta.donor}"
    conda "${moduleDir}/env.yaml"

    input:
    tuple val(meta), path(normal_bam), path(normal_bai), path(tumor_bam), path(tumor_bai)
    path reference_fasta
    path reference_index
    path reference_snps

    output:
    tuple val(meta), path("${meta.donor}_N.snpAF.tsv"), path("${meta.donor}_T.snpAF.tsv", optional: true), emit: snpAF

    script:

    def tumor_call = params.normal_tumor ? "alleleCounter -r ${reference_fasta} -l ${reference_snps} -b ${tumor_bam} -o ${meta.donor}_T.snpAF.tsv" : "touch ${meta.donor}_T.snpAF.tsv"

    """
    alleleCounter -r ${reference_fasta} -l ${reference_snps} -b ${normal_bam} -o ${meta.donor}_N.snpAF.tsv
    ${tumor_call}
    """
    stub:
    """
    touch ${meta.donor}_T.snpAF.tsv
    touch ${meta.donor}_N.snpAF.tsv
    """
}