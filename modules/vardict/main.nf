process VARDICT_CALL {

    publishDir "structural_varcalls/$meta.sample_name/cnkvit", mode: 'copy'
    tag "${meta.sample_name}"
    conda "${moduleDir}/env.yaml"

    input:
    tuple val(meta), path(bam), path(bai)
    path reference_fasta
    path reference_fasta_fai
    path organism_dna_panel


    output:
    tuple val(meta), path("results/vardict_SNV_${meta.sample_name}.vcf"), emit: vcfs

    script:

    """
    mkdir -p results
    vardict-java -G ${reference_fasta} -th $task.cpus -N ${meta.sample_name} -b ${bam} -c 1 -S 2 -E 3 -g 4 ${organism_dna_panel} | teststrandbias.R | var2vcf_valid.pl -m 7 -c 1 -N ${meta.sample_name} -f ${params.AF_threshold} > results/vardict_SNV_${meta.sample_name}.vcf
    sed -i '/TYPE=[DI][UNE][PVL]/d' results/vardict_SNV_${meta.sample_name}.vcf
    """

    stub:

    """
    mkdir -p results
    touch results/vardict_SNV_${meta.sample_name}.vcf
    """

}