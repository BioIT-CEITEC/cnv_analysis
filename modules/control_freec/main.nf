process CONTROL_FREEC {
    tag "${meta.donor}"
    publishDir "structural_varcalls/${meta.donor}/control_freec", mode: 'copy'
    conda "${moduleDir}/env.yaml"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bam_bai), path(normal_bam), path(normal_bam_bai)
    path reference_fasta
    path reference_index
    path snps_bed
    path(binned_bed)
    path(gc_profile)

    output:
    tuple val(meta), path("control_freec/*.CNV_varcalls.tsv"),  path("control_freec/config.txt"), emit: var_call

        script:

    def normal_config = params.normal_tumor ? """
    filedata = filedata.replace('#normal_X__', '')
    filedata = filedata.replace('X_input_control_bam_X', '${normal_bam}')
    """ : ""

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
    filedata = filedata.replace('X_input_bam_X', '${tumor_bam}')

    ${normal_config}

    with open('control_freec/config.txt', 'w') as file:
        file.write(filedata)
    CODE

    freec -conf control_freec/config.txt
    mv control_freec/${meta.donor}.bam_CNVs control_freec/${meta.donor}.CNV_varcalls.tsv
    """
}