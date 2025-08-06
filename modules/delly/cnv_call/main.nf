process CNV_CALLS_DELLY {

    publishDir "structural_varcalls/$meta.donor/delly", mode: 'copy'
    tag "Sample: ${meta.donor}"
    conda "${moduleDir}/../env.yaml"

    input:
    tuple val(meta), path(normal_bam), path(normal_bai), path(tumor_bam), path(tumor_bai)
    path reference_fasta
    path reference_fasta_fai
    path map_file // channel to the reference fasta fil
    tuple path(sv_calls), path(vcf), path(vcf_index)
    
  

    output:
    tuple val(meta), path("${meta.donor}.cnvs.bcf"), path("${meta.donor}.cnvs.bcf.csi"), emit: cnv_calls_delly
    
    script:
    
        """
        delly cnv \\
          -i 100000 -j 100000 -w 100000 \\
          -o ${meta.donor}.cnvs.bcf \\
          -g ${reference_fasta} \\
          -m ${map_file} \\
          -l ${sv_calls} \\
          ${normal_bam}
        """

    stub:
        """
        mkdir -p results
        touch results/target.bed
        touch results/antitarget.bed
        """
        
}