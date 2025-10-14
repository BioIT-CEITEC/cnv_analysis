process CONTROL_FREEC {
    tag "${meta.sample_name}"
    publishDir "structural_varcalls/", mode: 'copy'
    conda "${moduleDir}/env.yaml"

    input:
    tuple val(meta), path(bam), path(bam_bai)
    path reference_fasta
    path reference_index
    path snps_bed
    path(binned_bed)
    path(gc_profile)

    output:
    tuple val(meta), path("control_freec/*.CNV_varcalls.tsv"),  path("control_freec/config.txt"), path("control_freec/*_ratio.txt"), emit: var_call

    script:

    """
    mkdir -p control_freec
    python <<CODE
    import os

    template_file = "${moduleDir}/control_freec_config_template_WGS.txt"

    with open(template_file, 'r') as file:
        filedata = file.read()

    filedata = filedata.replace('X_window_X', str(${params.wgs_bin_size}))
    filedata = filedata.replace('X_ref_fai_X', '${reference_index}')
    filedata = filedata.replace('X_GC_profile_file_X', '${gc_profile}')
    filedata = filedata.replace('X_output_dir_X', 'control_freec/')
    filedata = filedata.replace('X_input_bam_X', '${bam}')

    with open('control_freec/config.txt', 'w') as file:
        file.write(filedata)
    CODE

    freec -conf control_freec/config.txt
    mv control_freec/${meta.sample_name}.bam_CNVs control_freec/${meta.sample_name}.CNV_varcalls.tsv
    """
}