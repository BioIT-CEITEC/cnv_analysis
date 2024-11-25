process COVERAGE_CALC {

    tag "${meta.id}"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    tuple path(binned_genome), path(gc_profile)
    path organism_dna_panel
    path organism_dict

    output:
    tuple val(meta), path("*.region_coverage.tsv"), emit: region_coverage

    script:

    def region_bed = params.tumor_normal_paired ? "${binned_genome}" : "${organism_dna_panel}"
    def normal_flag = params.tumor_normal_paired ? "bedtools coverage -sorted -a ${region_bed} -b ${normal_bam} -o ${organism_dict} > ${meta.id}.region_coverage.tsv" : ""

    """
    bedtools coverage -sorted -a ${region_bed} -b ${tumor_bam} -g ${organism_dict} > ${meta.id}.region_coverage.tsv
    ${normal_flag}
    """

    stub:
    """
    touch test.region_coverage.tsv
    """

    // TODO: include the threads to the process and the code
}