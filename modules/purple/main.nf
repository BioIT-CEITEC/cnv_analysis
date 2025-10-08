process PURPLE {
    tag "${meta.sample_name}"
    publishDir "structural_varcalls/${meta.sample_name}/purple", mode: 'copy'

    conda "${moduleDir}/env.yml"


    input:
    tuple val(meta), path(amber), path(cobalt), path(normal_vcf)
    path genome_fasta
    path genome_fai
    path genome_dict
    path gc_profile
    path germline_hotspots
    path driver_gene_panel
    path ensembl_data_resources
    path germline_del_freq

    output:
    tuple val(meta), path('purple/'), emit: purple_dir


    script:

    def genome_ver = params.assembly.replace("GRCh","")

    """

    purple \\
        -Xmx16G \\
        -reference ${meta.sample_name} \\
        -amber ${amber} \\
        -cobalt ${cobalt} \\
        -gc_profile ${gc_profile} \\
        -ref_genome_version ${genome_ver} \\
        -ref_genome ${genome_fasta} \\
        -ensembl_data_dir ${ensembl_data_resources}/ensembl_data/ \\
        -germline_hotspots ${germline_hotspots} \\
        -germline_del_freq_file ${germline_del_freq} \\
        -germline_vcf ${normal_vcf} \\
        -driver_gene_panel ${driver_gene_panel} \\
        -output_dir purple/

    """

    stub:
    """
    mkdir purple/
    touch purple/${meta.sample_name}.purple.cnv.gene.tsv
    touch purple/${meta.sample_name}.purple.cnv.somatic.tsv
    """
}
