process ALIGN_REGIONS {

    tag "${paired_bed.baseName}"
    publishDir "pseudogene_identification/aligned_regions", mode: 'copy'
    conda "${moduleDir}/../env.yaml"

    input:
    path reference
    path reference_fai
    path paired_bed

    output:
    path("diff_output"), emit: diff_tsvs

    script:
    """
    python ${projectDir}/bin/align_regions.py \
        --reference ${reference} \
        --bed ${paired_bed} \
        --output_dir diff_output/
    """
}
