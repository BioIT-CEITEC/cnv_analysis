process MERGE_REGIONS_DELLY {

    publishDir "structural_varcalls/all_samples/delly", mode: 'copy'
    tag "Merging prepared regions..."
    conda "${moduleDir}/../env.yaml"

    input:
    path prepared_regions // channel to the reference fasta fil

    output:
    path("dellySV.bcf"), emit: merged_regions_delly
    
    script:
    
        """
        delly merge \\
          -o dellySV.bcf \\
          ${prepared_regions}
        """

    stub:
        """
        mkdir -p results
        touch results/target.bed
        touch results/antitarget.bed
        """
        
}