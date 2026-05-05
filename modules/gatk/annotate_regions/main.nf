process ANNOTATE_REGIONS_GATK {

    conda "${moduleDir}/../env.yaml"

    input:
    path reference_fasta
    path reference_fasta_fai
    path lib_ROI
    path reference_dict
    path prepared_regions

    output:
    path("*.annotated.tsv"), emit: annotated_regions_gatk

    script:
    def panel = lib_ROI.toString().replace('.bed', '') 
    """
    gatk AnnotateIntervals \\
      -O ${panel}.annotated.tsv \\
      -R ${reference_fasta} \\
      -imr OVERLAPPING_ONLY \\
      -L ${prepared_regions} \\
    """

    stub:
    """
    touch placeholder.annotated.tsv
    """
}