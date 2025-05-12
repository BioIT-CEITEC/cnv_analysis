process VARDICT_CALL {

    publishDir "structural_varcalls/$meta.donor/cnkvit", mode: 'copy'
    tag "${meta.donor}"
    conda "${moduleDir}/env.yaml"

    input:
    tuple val(meta), path(normal_bam), path(normal_bai), path(tumor_bam), path(tumor_bai)
    path reference_fasta
    path reference_fasta_fai
    path organism_dna_panel


    output:
    tuple val(meta), path("results/vardict_SNV_${meta.donor}_N.vcf"), path("results/vardict_SNV_${meta.donor}_T.vcf", optional: true), emit: vcfs

    script:

    def tumor_cmd = params.normal_tumor ? "vardict-java -G ${reference_fasta} -th $task.cpus -N ${meta.donor}_T -b ${tumor_bam} -c 1 -S 2 -E 3 -g 4 ${organism_dna_panel} | teststrandbias.R | var2vcf_valid.pl -m 7 -c 1 -N ${meta.donor}_T -f ${params.AF_threshold} > results/vardict_SNV_${meta.donor}_T.vcf && sed -i '/TYPE=[DI][UNE][PVL]/d' results/vardict_SNV_${meta.donor}_T.vcf" : "touch results/vardict_SNV_${meta.donor}_T.vcf"

    """
    mkdir -p results
    vardict-java -G ${reference_fasta} -th $task.cpus -N ${meta.donor}_N -b ${normal_bam} -c 1 -S 2 -E 3 -g 4 ${organism_dna_panel} | teststrandbias.R | var2vcf_valid.pl -m 7 -c 1 -N ${meta.donor}_N -f ${params.AF_threshold} > results/vardict_SNV_${meta.donor}_N.vcf
    sed -i '/TYPE=[DI][UNE][PVL]/d' results/vardict_SNV_${meta.donor}_N.vcf
    ${tumor_cmd}
    """

    stub:

    """
    mkdir -p results
    touch results/vardict_SNV_${meta.donor}_N.vcf
    touch results/vardict_SNV_${meta.donor}_T.vcf
    """

}