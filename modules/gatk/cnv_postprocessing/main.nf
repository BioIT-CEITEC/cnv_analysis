process POSTPROCESSING_CNV_GATK {

    publishDir "structural_varcalls/${meta.sample_name}", mode: 'copy'
    tag "${meta.sample_name}"
    conda "${moduleDir}/../env.yaml"

    input:
    tuple val(meta), path(normal_bam), path(normal_bai), path(germline_calls)
    path cohort_data
    path ploidy_calls
    path reference_dict

    output:
    tuple val(meta), path("gatk/*_intervals_${meta.sample_name}.vcf.gz"),
    path("gatk/*_intervals_${meta.sample_name}.vcf.gz.tbi"),
    path("gatk/${meta.sample_name}_gatk.vcf.gz"),
    path("gatk/${meta.sample_name}_gatk.vcf.gz.tbi"),
    path("gatk/copy_ratios_${meta.sample_name}.tsv"), emit: postprocess_cnv_gatk

    script:
    def chrX = params.assembly ==~ /hg(19|38)/ ? "chrX" : "X"
    def chrY = params.assembly ==~ /hg(19|38)/ ? "chrY" : "Y"
    """
    mkdir -p gatk
    export XDG_CACHE_HOME=\$PWD/.cache
    mkdir -p \$XDG_CACHE_HOME/arviz
    gatk PostprocessGermlineCNVCalls \\
      --model-shard-path cohort-model \\
      --calls-shard-path germline-calls \\
      --allosomal-contig ${chrX} --allosomal-contig ${chrY} \\
      --contig-ploidy-calls ploidy-calls \\
      --sample-index 0 \\
      --output-genotyped-intervals gatk/genotyped_intervals_${meta.sample_name}.vcf.gz \\
      --output-genotyped-segments gatk/${meta.sample_name}_gatk.vcf.gz \\
      --output-denoised-copy-ratios gatk/copy_ratios_${meta.sample_name}.tsv \\
      --sequence-dictionary ${reference_dict}
    """

    stub:
    """
    mkdir -p gatk
    touch gatk/${meta.sample_name}.vcf.gz
    touch gatk/${meta.sample_name}.vcf.gz.tbi
    touch gatk/copy_ratios_${meta.sample_name}.tsv
    touch gatk/segments_intervals_${meta.sample_name}.vcf.gz
    touch gatk/segments_intervals_${meta.sample_name}.vcf.gz.tbi
    """
}