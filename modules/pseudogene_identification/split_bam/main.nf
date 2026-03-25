process SplitBam{
    tag "${params.bam}"
    conda './env.yml'
    publishDir path: "${params.output_dir}/${params.run_name}/${params.run_name}_split", mode: 'copy'
    input: 
    tuple path(tsv), path(bam)

    output: 
    tuple path("${bam.baseName}_gene.bam"), path("${bam.baseName}_pseudogene.bam"), path("${bam.baseName}_ambiguous.bam")

    script: 
    """

    awk -F '\t' '\$6 == "Gene" {OFS="\t"; print \$0}' "${tsv}" > gene_reads.tsv
    awk -F '\t' '\$6 == "Pseudogene" {OFS="\t"; print \$0}' "${tsv}" > pseudogene_reads.tsv
    awk -F '\t' '\$6 == "Ambiguous" {OFS="\t"; print \$0}' "${tsv}" > ambiguous_reads.tsv

    head gene_reads.tsv
    head pseudogene_reads.tsv
    head ambiguous_reads.tsv

    modify_bam.py --tsv gene_reads.tsv --bam ${bam} --output ${bam.baseName}_gene.bam
    modify_bam.py --tsv pseudogene_reads.tsv --bam ${bam} --output ${bam.baseName}_pseudogene.bam
    modify_bam.py --tsv ambiguous_reads.tsv --bam ${bam} --output ${bam.baseName}_ambiguous.bam

    """
}