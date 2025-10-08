process POSTPROCESSING_CNV_GATK {

    publishDir "structural_varcalls/$meta.donor/gatk", mode: 'copy'
    tag "sample: ${meta.sample_name}"
    conda "${moduleDir}/../env.yaml"

  input:

  tuple val(meta), path(normal_bam), path(normal_bai)
  path cohort_data // channel to the reference fasta file
  path ploidy_calls
  path reference_dict

    output:
    tuple val(meta), path("*_intervals_${meta.sample_name}.vcf.gz"),
    path("*_intervals_${meta.sample_name}.vcf.gz.tbi"),
    path("*_segments_${meta.sample_name}.vcf.gz"),
    path("*_segments_${meta.sample_name}.vcf.gz.tbi"),
    path("copy_ratios_${meta.sample_name}.tsv"), emit: postprocess_cnv_gatk

    script:

        """
        gatk PostprocessGermlineCNVCalls \\
          --model-shard-path cohort-model \\
          --calls-shard-path germline-calls \\
          --allosomal-contig X --allosomal-contig Y \\
          --contig-ploidy-calls ploidy-calls \\
          --sample-index 0 \\
          --output-genotyped-intervals genotyped_intervals_${meta.sample_name}.vcf.gz \\
          --output-genotyped-segments genotyped_segments_${meta.sample_name}.vcf.gz \\
          --output-denoised-copy-ratios copy_ratios_${meta.sample_name}.tsv \\
          --sequence-dictionary ${reference_dict}

        """

    stub:
        """
        mkdir -p results
        touch results/target.bed
        touch results/antitarget.bed
        """
}