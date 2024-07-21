process CONTROL_FREEC {

    input:
    tuple val(meta), path(bams)
    tuple path(binned_bed), path(gc_profile)

    output:
    tuple val(meta), path("CNV_varcalls.tsv"), path("control_freec/config.txt"), emit: var_call

    script:
    def prefix = "${meta.id}"
    
        if (!params.normal_tumor) {
        """
        [ ! -f  ${prefix}.bam ] && ln -s $bams ${prefix}_T.bam
        python <<CODE
    import os

    template_file = 'bin/control_freec_config_template_WGS.txt'

    with open(template_file, 'r') as file:
        filedata = file.read()

    filedata = filedata.replace('X_window_X', str('${params.window_size}'))
    filedata = filedata.replace('X_ref_fai_X', '${params.organism_fasta}.fai')
    filedata = filedata.replace('X_ref_X', 'GC_profile_${params.windows_size}.cnp')
    filedata = filedata.replace('X_output_dir_X', 'control_freec/')
    filedata = filedata.replace('X_input_bam_X', '${prefix}_T.bam)

    with open('control_freec/config.txt', 'w') as file:
        file.write(filedata)
    CODE

    freec -conf control_freec/config.txt
    mv ${prefix}.bam_CNVs CNV_varcalls.tsv
        """
    } else {
        """
        [ ! -f  ${prefix}_N.bam ] && ln -s ${bams[0]} ${prefix}_N.bam
        [ ! -f  ${prefix}_T.bam ] && ln -s ${bams[1]} ${prefix}_T.bam
        python <<CODE
    import os

    template_file = 'bin/control_freec_config_template_WGS.txt'

    with open(template_file, 'r') as file:
        filedata = file.read()

    filedata = filedata.replace('X_window_X', str(${params.window_size}))
    filedata = filedata.replace('X_ref_fai_X', '${params.organism_fasta}.fai')
    filedata = filedata.replace('X_ref_X', 'GC_profile_${params.windows_size}.cnp')
    filedata = filedata.replace('X_output_dir_X', 'control_freec/')
    filedata = filedata.replace('X_input_bam_X', '${prefix}_T.bam)
    filedata = filedata.replace('#normal_X__', '')
    filedata = filedata.replace('X_input_control_bam_X', ${prefix}_N.bam)

    with open('control_freec/config.txt', 'w') as file:
        file.write(filedata)
    CODE

    freec -conf control_freec/config.txt
    mv ${prefix}.bam_CNVs CNV_varcalls.tsv
        """
    }

    // TODO: include the threads to the process and the code
}


