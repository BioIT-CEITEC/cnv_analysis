process POSTPROCESSING_CNV_GATK {

    publishDir "structural_varcalls/$meta.donor/gatk", mode: 'copy'
    tag "sample: ${meta.donor}, index: ${meta.index_gatk}"
    conda "${moduleDir}/../env.yaml"

  input:

  tuple val(meta), path(normal_bam), path(normal_bai), path(tumor_bam), path(tumor_bai), path(germline_calls)
  path cohort_data // channel to the reference fasta file
  path ploidy_calls
  path reference_dict

    output:
    tuple val(meta), path("*_intervals_${meta.donor}.vcf.gz"),
    path("*_intervals_${meta.donor}.vcf.gz.tbi"),
    path("*_segments_${meta.donor}.vcf.gz"),
    path("*_segments_${meta.donor}.vcf.gz.tbi"),
    path("copy_ratios_${meta.donor}.tsv"), emit: postprocess_cnv_gatk

    script:

        """
        gatk PostprocessGermlineCNVCalls \\
          --model-shard-path cohort-model \\
          --calls-shard-path germline-calls \\
          --allosomal-contig X --allosomal-contig Y \\
          --contig-ploidy-calls ploidy-calls \\
          --sample-index ${meta.index_gatk} \\
          --output-genotyped-intervals genotyped_intervals_${meta.donor}.vcf.gz \\
          --output-genotyped-segments genotyped_segments_${meta.donor}.vcf.gz \\
          --output-denoised-copy-ratios copy_ratios_${meta.donor}.tsv \\
          --sequence-dictionary ${reference_dict}

        """

    stub:
        """
        mkdir -p results
        touch results/target.bed
        touch results/antitarget.bed
        """
}