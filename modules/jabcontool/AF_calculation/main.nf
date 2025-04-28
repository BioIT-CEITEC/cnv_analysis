process SNP_AF_CALC {

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path reference_fasta
    path reference_index
    path reference_snps

    output:
    tuple val(meta), path("${meta.tumor_id}.snpAF.tsv"), path("${meta.normal_id}*.snpAF.tsv", optional: true), emit: snpAF

    script:

    def prefix = "${meta.donor}"
    def normal_call = params.normal_tumor ? "alleleCounter -r ${reference_fasta} -l ${reference_snps} -b ${normal_bam} -o ${meta.normal_id}.snpAF.tsv" : ""

    """
    alleleCounter -r ${reference_fasta} -l ${reference_snps} -b ${tumor_bam} -o ${meta.tumor_id}.snpAF.tsv
    ${normal_call}
    """
    stub:
    """
    touch test_T.snpAF.tsv
    touch test_N.snpAF.tsv
    """
}