process COMBINE_BAM {
    tag "${meta.sample_name}_${gene_region}_${pseudogene_region}"
    publishDir "pseudogene/test", mode: 'copy'
    conda "${moduleDir}/../env.yaml"

    input:
    tuple val(meta), path(bam), path(bai), val(gene_region), path(gene_reference), val(pseudogene_region), path(pseudogene_reference)

    output:
    tuple val(meta), val(gene_region), val(pseudogene_region), path(gene_reference), path("${meta.sample_name}_combined_${gene_region}_${pseudogene_region}.bam"), emit: combined_bam

    script:
    """
    samtools view -b ${bam} ${gene_region} > gene_region_${gene_region}.bam
    samtools view -b ${bam} ${pseudogene_region} > pseudogene_region_${pseudogene_region}.bam

    samtools sort -o gene_region_sorted_${gene_region}.bam gene_region_${gene_region}.bam
    samtools sort -o pseudogene_region_sorted_${pseudogene_region}.bam pseudogene_region_${pseudogene_region}.bam

    samtools index gene_region_sorted_${gene_region}.bam
    samtools index pseudogene_region_sorted_${pseudogene_region}.bam

    samtools merge -f -o ${meta.sample_name}_combined_${gene_region}_${pseudogene_region}.bam gene_region_sorted_${gene_region}.bam pseudogene_region_sorted_${pseudogene_region}.bam
    samtools index ${meta.sample_name}_combined_${gene_region}_${pseudogene_region}.bam
    """
}