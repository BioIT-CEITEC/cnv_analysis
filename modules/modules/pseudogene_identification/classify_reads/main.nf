process CLASSIFY_READS {
    tag "${meta.sample_name}"

    conda "${moduleDir}/../env.yaml"

    input:
    tuple val(meta), path(realigned_dir), path(bam), path(bai)
    path diff_dir
    path region_bed

    output:
    tuple val(meta), path("${meta.sample_name}_classified_reads.tsv"), emit: classified_reads

    script:
    """
    python3 ${projectDir}/bin/classify_reads_by_base.py \
        --diff_dir      ${diff_dir} \
        --realigned_dir ${realigned_dir} \
        --bam_original  ${bam} \
        --bed           ${region_bed} \
        --sample        ${meta.sample_name} \
        --output        ${meta.sample_name}_classified_reads.tsv
    """

    stub:
    """
    touch ${meta.sample_name}_classified_reads.tsv
    """
}