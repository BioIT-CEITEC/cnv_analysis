process SNP_AF_CALC {

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path organism_fasta
    path organism_snps

    output:
    tuple val(meta), path("*_T.snpAF.tsv"), path("*_N.snpAF.tsv", optional: true), emit: snpAF

    script:

    def prefix = "${meta.id}"
    def normal_call = params.normal_tumor ? "alleleCounter -r ${organism_fasta} -l ${organism_snps} -b ${normal_bam} -o ${meta.id}_N.snpAF.tsv" : ""

    """
    alleleCounter -r $params.organism_fasta -l $params.organism_snps_tsv -b ${prefix}_N.bam -o ${meta.id}_T.snpAF.tsv
    ${normal_call}
    """
    stub:
    """
    touch test_T.snpAF.tsv
    touch test_N.snpAF.tsv
    """

    // TODO: include the threads to the process and the code
}