process Summarize {
    tag "${params.bam}"
    conda './env.yml'
    publishDir path: "${params.output_dir}/${params.run_name}", mode: 'copy'

    input: 
    tuple path(tsv), path(original_bam), path(params.gene_bed), path(params.pseudogene_bed)

    output:
    tuple path("${params.run_name}_report_all.csv"), path("${params.run_name}_plot.html"), path("${params.run_name}_unique_classification.tsv")

    script: 
    """
    summarize.py --classificaton_tsv ${tsv} --gene_bed ${params.gene_bed} --pseudogene_bed ${params.pseudogene_bed} --bam ${original_bam} --report_all ${params.run_name}_report_all.csv --plot_output ${params.run_name}_plot.html --output_classification ${params.run_name}_unique_classification.tsv

    """
}