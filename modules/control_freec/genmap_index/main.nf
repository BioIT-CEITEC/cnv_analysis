process GENMAP_INDEX_FREEC {
    conda "${moduleDir}/../env.yaml"

    input:
    path reference_fasta
    path reference_fasta_fai

    output:
    path("genmap_index/*"), emit: genmap_index
    path("freec_chrLen.txt"), emit: freec_chrlen

    script:
    """
    genmap index -F ${reference_fasta} -I  genmap_index/
    cut -f1,2 ${reference_fasta_fai} > freec_chrLen.txt
    """

    stub:
    """
    touch genmap_index/placeholder.txt
    touch freec_chrLen.txt
    """
}