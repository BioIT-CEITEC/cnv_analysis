process ForceMap {
    tag "${gene_region}_${pseudogene_region}"
    conda './env.yml'
    input:
    tuple val(gene_region), val(pseudogene_region), path(fasta), path(bam)

    output:
    tuple val(gene_region), val(pseudogene_region), path("${bam.baseName}_force_mapped.bam")

    script:
    """
    samtools fastq -n ${bam} > ${bam.baseName}_reads.fq
    samtools faidx ${fasta}
    bwa index ${fasta}
    bwa mem -T 0 -A 1 -B 1 -O 1 -E 1 -t 4 ${fasta} ${bam.baseName}_reads.fq | samtools view -b - | samtools sort -o ${bam.baseName}_tmp.bam && samtools index ${bam.baseName}_tmp.bam

    reindex_remapped_bam.py --bam ${bam.baseName}_tmp.bam --output ${bam.baseName}_force_mapped.bam --header ${bam}
    samtools index ${bam.baseName}_force_mapped.bam
    """
}