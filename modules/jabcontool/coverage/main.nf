process COVERAGE_CALC {

    tag "${meta.donor}"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path organism_reference
    path reference_index

    output:
    tuple val(meta), path("${meta.tumor_id}.region_coverage.tsv"), path("${meta.normal_id}.region.coverage.tsv", optional: true), emit: region_coverage

    script:

    def normal_flag = params.tumor_normal_paired ? "bedtools coverage -sorted -a ${organism_reference} -b ${normal_bam} -o ${reference_index} > ${meta.normal_id}.region_coverage.tsv" : ""

    """
    bedtools coverage -sorted -a ${organism_reference} -b ${tumor_bam} -g ${reference_index} > ${meta.tumor_id}.region_coverage.tsv
    ${normal_flag}
    """

    stub:
    """
    touch test.region_coverage.tsv
    """
}