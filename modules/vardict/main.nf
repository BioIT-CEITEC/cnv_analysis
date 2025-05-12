process VARDICT_CALL {

    publishDir "structural_varcalls/${meta.donor}/cnkvit", mode: 'copy'
    tag "${meta.donor}"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path reference_fasta
    path reference_fasta_fai
    path organism_dna_panel


    output:
    tuple val(meta), path("*.${meta.tumor_id}.vcf"), path("*.${meta.normal_id}.vcf", optional: true), emit: vcfs

    script:

    def normal_cmd = params.normal_tumor ? "vardict -java -G ${reference_fasta} -th $task.cpus -N ${meta.normal_id} -b ${normal_bam} -c 1 -S 2 -E 3 -g 4 ${organism_dna_panel} | teststrandbias.R | var2vcf_valid.pl -m 7 -c 1 -N ${meta.normal_id} -f ${params.AF_threshold} > results/vardict_SNV_${meta.normal_id}.vcf && sed -i '/TYPE=[DI][UNE][PVL]/d' results/vardict_SNV_${meta.normal_id}.vcf" : ""

    """
    mkdir -p results
    vardict -java -G ${reference_fasta} -th $task.cpus -N ${meta.tumor_id} -b ${tumor_bam} -c 1 -S 2 -E 3 -g 4 ${organism_dna_panel} | teststrandbias.R | var2vcf_valid.pl -m 7 -c 1 -N ${meta.tumor_id} -f ${params.AF_threshold} > results/vardict_SNV_${meta.tumor_id}.vcf
    sed -i '/TYPE=[DI][UNE][PVL]/d' results/vardict_SNV_${meta.tumor_id}.vcf
    ${normal_cmd}
    """

    stub:

    """
    mkdir -p results
    touch results/vardict_SNV_test_T.vcf
    touch results/vardict_SNV_test_N.vcf
    """

}