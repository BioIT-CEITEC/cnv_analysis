process AMBER {
    tag "${meta.sample_name}"

    conda "${moduleDir}/env.yml"

    input:
    tuple val(meta), path(normal_bam), path(normal_bai)
    path heterozygous_sites
    path target_region_bed

    output:
    tuple val(meta), path('amber/'), emit: amber_dir

    script:

    def target_regions_bed_arg = params.lib_ROI != "wgs" ? "-target_regions_bed ${target_region_bed}" : ''
    def genome_ver = params.assembly.replace("GRCh", "")

    """
    amber \\
        -Xmx16G \\
        -reference ${meta.sample_name}_N \\
        -reference_bam ${normal_bam} \\
        ${target_regions_bed_arg} \\
        -ref_genome_version ${genome_ver} \\
        -loci ${heterozygous_sites} \\
        -output_dir amber/

    """

    stub:
    """
    mkdir -p amber/
    touch amber/placeholder
    """
}
