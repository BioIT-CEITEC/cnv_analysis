process READ_PROCESSING {
    tag "${meta.sample_name}_${gene_region}_${pseudogene_region}"
    publishDir "pseudogene/test", mode: 'copy'
    conda "${moduleDir}/../env.yaml"

    input:
    tuple val(meta), path(bam), val(gene_region), val(pseudogene_region), path(tsv)

    output:
    tuple val(meta), val(gene_region), val(pseudogene_region), path("${meta.sample_name}_${gene_region}_${pseudogene_region}_reads.tsv"), emit: processed_reads

    //${params.min_mapq} remember to add this
    script:
    """
    samtools index ${bam}
    python ${projectDir}/bin/reads_processing.py --tsv ${tsv} --bam ${bam} --output ${meta.sample_name}_${gene_region}_${pseudogene_region}_reads.tsv --min_mapq 10 
    """
}