process COBALT {
    tag "${meta.donor}"

    conda "${moduleDir}/env.yml"


    input:
    tuple val(meta), path(normal_bam), path(normal_bai), path(tumor_bam), path(tumor_bai)
    path gc_profile
    path diploid_regions
    path organism_panel

    output:
    tuple val(meta), path('cobalt/'), emit: cobalt_dir

    script:


def target_region_arg = params.lib_ROI != "wgs" ? "-target_regions_bed ${organism_panel}" : ''

    """
    cobalt \\
        -Xmx16G \\
        -reference ${meta.donor}_N \\
        -reference_bam ${normal_bam} \\
        -gc_profile ${gc_profile} \\
        -output_dir cobalt/

    """

    stub:
    """
    mkdir -p cobalt/
    touch cobalt/placeholder

    """
}