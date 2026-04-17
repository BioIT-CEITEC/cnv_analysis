process EXTRACT_REGIONS {
    tag "${meta.sample_id}"
    conda "${moduleDir}/../env.yaml"

    input:
    tuple val(meta), path(bam), path(bai)
    path region_bed

    output:
    tuple val(meta), path("${meta.sample_id}_extracted_reads/"), path(bam), path(bai), emit: extracted_reads

    script:
    """
    python3 ${projectDir}/bin/extract_regions.py --bam ${bam} --bed ${region_bed} --outdir ${meta.sample_id}_extracted_reads/
    """
}