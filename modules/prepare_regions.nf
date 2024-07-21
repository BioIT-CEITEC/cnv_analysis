process PREPARE_REGIONS_CNVKIT {

    publishDir "structural_varcalls/all_samples/cnvkit", mode: 'copy'

    conda (params.conda_enabled ? "bioconda::bioconductor-copynumber" : null)

    input:
    ch_reference_fasta

    output:
    tuple path("*.target.bed"), path("*.antitarget.bed"), emit: cnvkit_beds

    script:

        """
        cnvkit.py access $params.organism_fasta -o reference_bed.bed
        cnvkit.py autobin $bams -t $regions -g reference_bed.bed
        mv ${params.lib_ROI}.target.bed target.bed
        mv ${params.lib_ROI}.antitarget.bed antitarget.bed
        """
}

process GET_COVERAGE_CNVKIT {

    tag "$meta.id"
    publishDir "structural_varcalls/$meta.id/cnvkit", mode: 'copy'

    input:
    tuple val(meta), path(bams), path(target), path(antitarget)

    output:
    tuple val(meta), path("normal.targetcoverage.cnn"), path("normal.antitargetcoverage.cnn"), path("tumor.targetcoverage.cnn"), path("tumor.antitargetcoverage.cnn"), emit: coverage

    script:

    def prefix = "${meta.id}"


        if (!params.normal_tumor) {
        """
        [ ! -f  ${prefix}.bam ] && ln -s $bams ${prefix}_T.bam
        cnvkit.py coverage ${prefix}_T.bam $target -o tumor.targetcoverage.cnn
        cnvkit.py coverage ${prefix}_T.bam $antitarget -o tumor.antitargetcoverage.cnn
        touch normal.targetcoverage.cnn
        touch normal.antitargetcoverage.cnn
        """
    } else {
        """
        [ ! -f  ${prefix}_N.bam ] && ln -s ${bams[0]} ${prefix}_N.bam
        [ ! -f  ${prefix}_T.bam ] && ln -s ${bams[1]} ${prefix}_T.bam
        cnvkit.py coverage ${prefix}_N.bam $target -o normal.targetcoverage.cnn
        cnvkit.py coverage ${prefix}_N.bam $antitarget -o normal.antitargetcoverage.cnn
        cnvkit.py coverage ${prefix}_T.bam $target -o tumor.targetcoverage.cnn
        cnvkit.py coverage ${prefix}_T.bam $antitarget -o tumor.antitargetcoverage.cnn
        """
    }

}

process REFERENCE_CNVKIT {

    publishDir "structural_varcalls/reference/cnvkit", mode: 'copy'

    input:
    path(reference)
    path(normal_cov)
    path(antitarget_cov)

    output:
    path("*.cnn"), emit: reference

    script:

        if (!params.normal_tumor) {
        """
        cnvkit.py reference --fasta $params.organism_fasta \
            -o normal_reference.cnn \
            -t $target_cov \
            -a $antitarget_cov
        """
    } else {
        """
        cnvkit.py reference $normal \
            --fasta $params.organism_fasta \
            -o normal_reference.cnn
        """
    }

}

process CNVKIT_CALL {

    publishDir "structural_varcalls/$meta.id/cnvkit", mode: 'copy'

    input:
    tuple val(meta), path(normal_target_cov), path(normal_antitarget_cov), path(tumor_target_cov), path(tumor_antitarget_cov), path(vcfs)
    path(reference)

    output:
    tuple val(meta), path("*.cnr"), path("*_cov.cns"), emit: fix_seg
    tuple val(meta), path("CNV_calls.cns"), path("cnvkit_diagram.pdf"), path("cnvkit_scatter.png"), emit: cnvkit

    script:

    def prefix = "${meta.id}"

    """
    cnvkit.py fix $tumor_target_cov $tumor_antitarget_cov $reference -o fixed_cov.cnr
    cnvkit.py segment fixed_cov.cnr -o segmented_cov.cns
    cnvkit.py call -y -m clonal segmented_cov.cns -o CNV_calls.cns --purity 0.5
    cnvkit.py diagram fixed_cov.cnr -s CNV_calls.cns -o cnvkit_diagram.pdf
    cnvkit.py scatter fixed_cov.cnr -s CNV_calls.cns -v vardict_SNV_${prefix}_T.vcf -o cnvkit_scatter.png
    """
}

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