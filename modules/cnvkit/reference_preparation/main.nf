process REFERENCE_CNVKIT {

    publishDir "structural_varcalls/reference/cnvkit", mode: 'copy'

    input:
    path reference
    path normal_cov
    path antitarget_cov

    output:
    path("*.cnn"), emit: reference

    script:

        if !(params.normal_tumor) {
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