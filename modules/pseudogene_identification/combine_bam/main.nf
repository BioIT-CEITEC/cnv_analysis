process CombineBam {
    tag "${gene_region}_${pseudogene_region}"
    conda './env.yml'
    input:
    tuple val(gene_region), val(pseudogene_region), file(bam)

    output:
    tuple val(gene_region), val(pseudogene_region), path("combined_${gene_region}_${pseudogene_region}.bam")

    script:
    """
    samtools view -b ${bam} ${gene_region} > gene_region_${gene_region}.bam
    samtools view -b ${bam} ${pseudogene_region} > pseudogene_region_${pseudogene_region}.bam

    samtools sort -o gene_region_sorted_${gene_region}.bam gene_region_${gene_region}.bam
    samtools sort -o pseudogene_region_sorted_${pseudogene_region}.bam pseudogene_region_${pseudogene_region}.bam

    samtools index gene_region_sorted_${gene_region}.bam
    samtools index pseudogene_region_sorted_${pseudogene_region}.bam

    samtools merge -f -o combined_${gene_region}_${pseudogene_region}.bam gene_region_sorted_${gene_region}.bam pseudogene_region_sorted_${pseudogene_region}.bam
    samtools index combined_${gene_region}_${pseudogene_region}.bam
    """
}