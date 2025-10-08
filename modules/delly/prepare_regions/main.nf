process PREPARE_REGIONS_DELLY {

    publishDir "structural_varcalls/$meta.sample_name/delly", mode: 'copy'
    tag "Sample: ${meta.sample_name}"
    conda "${moduleDir}/../env.yaml"

    input:
    tuple val(meta), path(normal_bam), path(normal_bai)
    path reference_fasta // channel to the reference fasta fil
    path reference_fasta_fai
    path reference_exclusion

    output:
    tuple val(meta), path("${meta.sample_name}.bcf"), emit: prepared_regions_delly

    script:

        """
        delly call \\
          -x ${reference_exclusion} \\
          -o ${meta.sample_name}.bcf \\
          -g ${reference_fasta} \\
          ${normal_bam}
        """

    stub:
        """
        mkdir -p results
        touch results/target.bed
        touch results/antitarget.bed
        """
}