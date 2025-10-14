process COVERAGE_CALC {

    tag "${meta.sample_name}"
    conda "${moduleDir}/env.yaml"

    input:
    tuple val(meta), path(bam), path(bam_bai)
    path organism_reference
    path reference_index

    output:
    tuple val(meta), path("${meta.sample_name}.region_coverage.tsv"), emit: region_coverage

    script:

    //def hasNormals = params.panel_of_normals ?: false
    //def tumor_flag = hasNormals ? "bedtools coverage -sorted -a ${organism_reference} -b ${tumor_bam} -o ${reference_index} > ${meta.sample_name}_T.region_coverage.tsv" : "touch ${meta.sample_name}_T.region.coverage.tsv"

    """
    bedtools coverage -sorted -a ${organism_reference} -b ${bam} -g ${reference_index} > ${meta.sample_name}.region_coverage.tsv
    """

    stub:
    """
    touch ${meta.sample_name}.region_coverage.tsv
    """
}