process VARDICT_CALL {

    publishDir "structural_varcalls/$meta.id/cnkvit", mode: 'copy'
    tag "$meta.id"

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("*.vcf"), emit: vcfs

    script:

    def prefix = "${meta.id}"

        if (!params.normal_tumor) {
        """
        [ ! -f  ${prefix}.bam ] && ln -s $bams ${prefix}_T.bam
        vardict -java -G $params.organism_fasta -th $task.cpus -N ${prefix}_T -b ${prefix}_T.bam -c 1 -S 2 -E 3 -g 4 $params.organism_dna_panel | teststrandbias.R | var2vcf_valid.pl -m 7 -c 1 -N ${prefix}_T -f $params.AF_threshold > vardict_SNV_${prefix}_T.vcf
        sed -i '/TYPE=[DI][UNE][PVL]/d' vardict_SNV_${prefix}_T.vcf
        """
    } else {
        """
        [ ! -f  ${prefix}_N.bam ] && ln -s ${bams[0]} ${prefix}_N.bam
        [ ! -f  ${prefix}_T.bam ] && ln -s ${bams[1]} ${prefix}_T.bam
        vardict -java -G $params.organism_fasta -th $task.cpus -N ${prefix}_N -b ${prefix}_N.bam -c 1 -S 2 -E 3 -g 4 $params.organism_dna_panel | teststrandbias.R | var2vcf_valid.pl -m 7 -c 1 -N ${prefix}_N -f $params.AF_threshold > vardict_SNV_${prefix}_N.vcf
        vardict -java -G $params.organism_fasta -th $task.cpus -N ${prefix}_T -b ${prefix}_T.bam -c 1 -S 2 -E 3 -g 4 $params.organism_dna_panel | teststrandbias.R | var2vcf_valid.pl -m 7 -c 1 -N ${prefix}_T -f $params.AF_threshold > vardict_SNV_${prefix}_T.vcf

        sed -i '/TYPE=[DI][UNE][PVL]/d' vardict_SNV_${prefix}_N.vcf
        sed -i '/TYPE=[DI][UNE][PVL]/d' vardict_SNV_${prefix}_T.vcf
        """
    }
}