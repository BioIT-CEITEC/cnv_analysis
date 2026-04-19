process REALIGN_READS {
    tag "${meta.sample_name}"
    conda "${moduleDir}/../env.yaml"

    input:
    tuple val(meta), path(extracted_dir), path(bam), path(bai)
    path reference_fasta
    path reference_fasta_fai
    path region_bed

    output:
    tuple val(meta), path("${meta.sample_name}_realigned_reads/"), path(bam), path(bai), emit: realigned_reads

    script:
    """
    python3 ${projectDir}/bin/realign_specific.py --bam_dir ${extracted_dir} --ref ${reference_fasta} --bed ${region_bed} --outdir ${meta.sample_name}_realigned_reads --flank 300
    """
}