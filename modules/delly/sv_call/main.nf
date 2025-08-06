process SV_CALLS_DELLY {

    publishDir "structural_varcalls/$meta.donor/delly", mode: 'copy'
    tag "Sample: ${meta.donor}"
    conda "${moduleDir}/../env.yaml"

    input:
    tuple val(meta), path(normal_bam), path(normal_bai), path(tumor_bam), path(tumor_bai)
    path reference_fasta
    path reference_fasta_fai
    path merged_regions // channel to the reference fasta fil
    path excluded_regions

    output:
    tuple val(meta), path("${meta.donor}.genotypeSV.bcf"), path("${meta.donor}.genotypeSV.bcf.csi"), emit: sv_genotype_delly
    
    script:
    
        """
        delly call \\
          -g ${reference_fasta} \\
          -v ${merged_regions} \\
          -o ${meta.donor}.genotypeSV.bcf \\
          -x ${excluded_regions} \\
          ${normal_bam}
        """

    stub:
        """
        mkdir -p results
        touch results/target.bed
        touch results/antitarget.bed
        """
        
}