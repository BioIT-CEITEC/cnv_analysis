process CLASSIFY_READS {
    tag "${meta.sample_id}"
    publishDir "pseudogene_identification/", mode: 'copy'
    conda "${moduleDir}/../env.yaml"

    input:
    tuple val(meta), path(realigned_dir), path(bam), path(bai)
    path region_bed

    output:
    tuple val(meta), path("${meta.sample_id}_classified_reads.tsv"), emit: classified_reads

    script:
    """
    python3 ${projectDir}/bin/reads_classification.py --bam_original ${bam} --realigned_dir ${realigned_dir} --bed ${region_bed} --out  ${meta.sample_id}_classified_reads.tsv
    """
}