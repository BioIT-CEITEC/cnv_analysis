process PREPARE_REGIONS_DELLY {

    publishDir "structural_varcalls/$meta.donor/delly", mode: 'copy'
    tag "Sample: ${meta.donor}"
    conda "${moduleDir}/../env.yaml"

    input:
    tuple val(meta), path(normal_bam), path(normal_bai), path(tumor_bam), path(tumor_bai)
    path reference_fasta // channel to the reference fasta fil
    path reference_fasta_fai
    path reference_exclusion
    
    output:
    tuple val(meta), path("${meta.donor}.bcf"), emit: prepared_regions_delly
    
    script:
    
        """
        delly call \\
          -x ${reference_exclusion} \\
          -o ${meta.donor}.bcf \\
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