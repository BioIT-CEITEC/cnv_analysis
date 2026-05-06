process PREPARE_REGIONS_GATK {

    conda "${moduleDir}/../env.yaml"

    input:
    path reference_fasta
    path reference_fasta_fai
    path lib_ROI
    path reference_dict

    output:
    path("*.interval_list"), emit: prepared_regions_gatk

    script:
    def panel = lib_ROI.toString().replace('.bed', '') 
    def flags = panel != "wgs" ? "-L ${lib_ROI} --bin-length 0" : "--padding 0"  
    """
    gatk PreprocessIntervals \\
      -O ${panel}.interval_list \\
      -R ${reference_fasta} \\
      -imr OVERLAPPING_ONLY \\
      ${flags} \\
    """

    stub:
    """
    touch placeholder.interval_list
    """
}