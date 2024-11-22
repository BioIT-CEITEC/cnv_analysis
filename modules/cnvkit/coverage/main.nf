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