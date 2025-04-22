process VARDICT_CALL {

    publishDir "structural_varcalls/$meta.id/cnkvit", mode: 'copy'
    tag "${meta.id}"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path reference_fasta
    path organism_dna_panel


    output:
    tuple val(meta), path("*_T.vcf"), path("*_N.vcf", optional: true), emit: vcfs

    script:

    def normal_flag = params.normal_tumor ? "" : " -N ${meta.id}_N"

    """
    mkdir -p results
    vardict -java -G ${reference_fasta} -th $task.cpus -N ${meta.id}_T -b ${tumor_bam} -c 1 -S 2 -E 3 -g 4 ${organism_dna_panel} | teststrandbias.R | var2vcf_valid.pl -m 7 -c 1 -N ${meta.id}_T -f ${params.AF_threshold} > results/vardict_SNV_${meta.id}_T.vcf
    sed -i '/TYPE=[DI][UNE][PVL]/d' results/vardict_SNV_${meta.id}_T.vcf
    """

    if (!params.normal_tumor) {
    """
    vardict -java -G ${reference_fasta} -th $task.cpus -N ${meta.id}_N -b ${normal_bam} -c 1 -S 2 -E 3 -g 4 ${organism_dna_panel} | teststrandbias.R | var2vcf_valid.pl -m 7 -c 1 -N ${meta.id}_N -f ${params.AF_threshold} > results/vardict_SNV_${meta.id}_N.vcf
    sed -i '/TYPE=[DI][UNE][PVL]/d' results/vardict_SNV_${meta.id}_N.vcf
    """
    }

    stub:

    """
    mkdir -p results
    touch results/vardict_SNV_test_T.vcf
    touch results/vardict_SNV_test_N.vcf
    """

}