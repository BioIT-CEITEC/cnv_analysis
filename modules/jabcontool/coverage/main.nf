process COVERAGE_CALC {

    tag "${meta.donor}"
    conda "${moduleDir}/env.yaml"

    input:
    tuple val(meta), path(normal_bam), path(normal_bai), path(tumor_bam), path(tumor_bai)
    path organism_reference
    path reference_index

    output:
    tuple val(meta), path("${meta.donor}_N.region_coverage.tsv"), path("${meta.donor}_T.region.coverage.tsv", optional: true), emit: region_coverage

    script:

    def tumor_flag = params.normal_tumor ? "bedtools coverage -sorted -a ${organism_reference} -b ${tumor_bam} -o ${reference_index} > ${meta.donor}_T.region_coverage.tsv" : "touch ${meta.donor}_T.region.coverage.tsv"

    """
    bedtools coverage -sorted -a ${organism_reference} -b ${normal_bam} -g ${reference_index} > ${meta.donor}_N.region_coverage.tsv
    ${tumor_flag}
    """

    stub:
    """
    touch ${meta.donor}_N.region_coverage.tsv
    touch ${meta.donor}_T.region_coverage.tsv
    """
}