process ALIGN_REGIONS {

    tag "${paired_bed.baseName}"
    publishDir "pseudogene_analysis_results", mode: 'copy'
    conda "${moduleDir}/../env.yaml"

    input:
    path reference
    path reference_fai
    path paired_bed

    output:
    path("aligned_regions"), emit: diff_tsvs

    script:
    """
    python ${projectDir}/bin/align_regions.py \
        --reference ${reference} \
        --bed ${paired_bed} \
        --output_dir aligned_regions/
    """

    stub:
    """
    mkdir -p aligned_regions
    """
}
