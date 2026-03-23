process POSTPROCESSING_CNV_GATK {

    publishDir "structural_varcalls/gatk", mode: 'copy'
    tag "sample: ${meta.sample_name}"
    conda "${moduleDir}/../env.yaml"

  input:

  tuple val(meta), path(normal_bam), path(normal_bai), path(germline_calls)
  path cohort_data // channel to the reference fasta file
  path ploidy_calls
  path reference_dict

    output:
    tuple val(meta), path("calls/*_intervals_${meta.sample_name}.vcf.gz"),
    path("calls/*_intervals_${meta.sample_name}.vcf.gz.tbi"),
    path("calls/${meta.sample_name}_gatk.vcf.gz"),
    path("calls/${meta.sample_name}_gatk.vcf.gz.tbi"),
    path("calls/copy_ratios_${meta.sample_name}.tsv"), emit: postprocess_cnv_gatk

    script:

        """
        mkdir -p calls
        gatk PostprocessGermlineCNVCalls \\
          --model-shard-path cohort-model \\
          --calls-shard-path germline-calls \\
          --allosomal-contig X --allosomal-contig Y \\
          --contig-ploidy-calls ploidy-calls \\
          --sample-index 0 \\
          --output-genotyped-intervals calls/genotyped_intervals_${meta.sample_name}.vcf.gz \\
          --output-genotyped-segments calls/${meta.sample_name}_gatk.vcf.gz \\
          --output-denoised-copy-ratios calls/copy_ratios_${meta.sample_name}.tsv \\
          --sequence-dictionary ${reference_dict}

        """

    stub:
        """
        mkdir -p results
        touch results/target.bed
        touch results/antitarget.bed
        """
}