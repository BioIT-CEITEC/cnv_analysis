process FORCEMAP {
    tag "${meta.sample_name}_${gene_region}_${pseudogene_region}"
    publishDir "pseudogene/test", mode: 'copy'
    conda "${moduleDir}/../env.yaml"

    input:
    tuple val(meta), val(gene_region), val(pseudogene_region), path(fasta), path(bam)

    output:
    tuple val(meta), val(gene_region), val(pseudogene_region), path("${meta.sample_name}_force_mapped.bam"), emit: force_mapped_bam

    script:
    """
    samtools fastq -n ${bam} > ${meta.sample_name}_reads.fq
    samtools faidx ${fasta}
    bwa index ${fasta}
    bwa mem -T 0 -A 1 -B 1 -O 1 -E 1 -t 4 ${fasta} ${meta.sample_name}_reads.fq | samtools view -b - | samtools sort -o ${meta.sample_name}_tmp.bam && samtools index ${meta.sample_name}_tmp.bam

    python ${projectDir}/bin/reindex_remapped_bam.py --bam ${meta.sample_name}_tmp.bam --output ${meta.sample_name}_force_mapped.bam --header ${bam}
    samtools index ${meta.sample_name}_force_mapped.bam
    """
}