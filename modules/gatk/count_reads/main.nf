process COUNT_READS_GATK {
    tag "${meta.sample_name}"
    conda "${moduleDir}/../env.yaml"

    input:
    tuple val(meta), path(normal_bam), path(normal_bai)
    path reference_fasta // channel to the reference fasta file
    path reference_fasta_fai
    path lib_ROI // channel to the regions of interest bed file
    path reference_dict
    path interval_list

    output:
    tuple val(meta), path("*.tsv"), emit: read_counts_gatk

    script:


        """
      gatk CollectReadCounts \
              -L ${interval_list} \
              -R ${reference_fasta} \
              -imr OVERLAPPING_ONLY \
              -I ${normal_bam} \
              --format TSV \
              -O ${meta.sample_name}.tsv

        """

    stub:
        """
        mkdir -p results
        touch results/target.bed
        touch results/antitarget.bed
        """

}
