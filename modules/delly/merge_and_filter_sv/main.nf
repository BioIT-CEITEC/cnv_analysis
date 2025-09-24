process MERGE_AND_FILTER_SV_DELLY {

    publishDir "structural_varcalls/all_samples/delly", mode: 'copy'
    tag "Merging and filtering SV calls"
    conda "${moduleDir}/../env.yaml"

    input:
    path genotype_calls
    path genotype_index

    output:
    tuple path("merged_genotypes.bcf"), path("germline_sv_calls.vcf.gz"), path("germline_sv_calls.vcf.gz.csi"), emit: final_sv_delly

    script:

        """
        bcftools merge \\
          -m id \\
          -O b \\
          -o merged_genotypes.bcf \\
          ${genotype_calls}

        bcftools index merged_genotypes.bcf

        delly filter \\
          -f germline \\
          -o germline_sv_calls.bcf \\
          merged_genotypes.bcf

        bcftools view germline_sv_calls.bcf -O z -o germline_sv_calls.vcf.gz
        bcftools index germline_sv_calls.vcf.gz
        """

    stub:
        """
        mkdir -p results
        touch results/target.bed
        touch results/antitarget.bed
        """
}