process READ_CLASSIFICATION {
    tag "${gene_region}_${pseudogene_region}"
    publishDir "pseudogene/test", mode: 'copy'
    conda "${moduleDir}/../env.yaml"

    input:
    tuple val(meta), val(gene_region), val(pseudogene_region), path(tsv)

    output:
    tuple val(meta),val(gene_region), val(pseudogene_region), path("${meta.sample_name}_${gene_region}_${pseudogene_region}_classified.tsv"), path("${meta.sample_name}_${gene_region}_${pseudogene_region}_depth.tsv"), path("${meta.sample_name}_${gene_region}_${pseudogene_region}_report.csv"), path("${meta.sample_name}_${gene_region}_${pseudogene_region}_coverage_plot.png"), path("*.png")

    script:
    //${params.threshold} remember to add this
    // metric is "threshold"
    """
    python ${projectDir}/bin/classify_reads.py --tsv ${tsv} --output "${meta.sample_name}_${gene_region}_${pseudogene_region}_classified.tsv" --threshold 0.7 --metric threshold
    python ${projectDir}/bin/calculate_depth_region.py --class_tsv "${meta.sample_name}_${gene_region}_${pseudogene_region}_classified.tsv" --positions ${tsv} --output_coverage "${meta.sample_name}_${gene_region}_${pseudogene_region}_depth.tsv"  --report "${meta.sample_name}_${gene_region}_${pseudogene_region}_report.csv" --output_plot ${meta.sample_name}_${gene_region}_${pseudogene_region}_coverage_plot.png
    """
} 