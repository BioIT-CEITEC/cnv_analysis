process ReadClassification {
    tag "${gene_region}_${pseudogene_region}"
    conda './env.yml'
    publishDir path: "${params.output_dir}/${params.run_name}/${gene_region}_${pseudogene_region}", mode: 'copy'

    input:
    tuple val(gene_region), val(pseudogene_region), path(tsv)

    output:
    tuple val(gene_region), val(pseudogene_region), path("${tsv.baseName}_classified.tsv"), path("${tsv.baseName}_depth.tsv"), path("${tsv.baseName}_report.csv"), path("${tsv.baseName}_coverage_plot.png"), path("*.png")

    script:
    """
    classify_reads.py --tsv ${tsv} --output "${tsv.baseName}_classified.tsv" --threshold ${params.threshold} --metric ${params.metric}
    calculate_depth_region.py --class_tsv "${tsv.baseName}_classified.tsv" --positions ${tsv} --output_coverage "${tsv.baseName}_depth.tsv"  --report "${tsv.baseName}_report.csv" --output_plot ${tsv.baseName}_coverage_plot.png
    """
} 