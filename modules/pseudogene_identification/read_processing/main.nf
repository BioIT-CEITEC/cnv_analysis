process ReadProcessing {
    tag "${gene_region}_${pseudogene_region}"
    conda './env.yml'
    input:
    tuple val(gene_region), val(pseudogene_region), path(tsv), path(bam)

    output:
    tuple val(gene_region), val(pseudogene_region), path("${bam.baseName}_reads.tsv")

    script:
    """
    samtools index ${bam}
    reads_processing.py --tsv ${tsv} --bam ${bam} --output ${bam.baseName}_reads.tsv --min_mapq ${params.min_mapq}
    """
}