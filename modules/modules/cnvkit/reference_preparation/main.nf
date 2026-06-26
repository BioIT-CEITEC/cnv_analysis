process REFERENCE_CNVKIT {

    conda "${moduleDir}/../env.yaml"

    input:
    path reference_fasta
    path normal_coverage
    path sample_coverage
    tuple path(target), path(antitarget)
    val sample_number

    output:
    path("reference/normal_reference.cnn"), emit: cnvkit_reference

    script:



    if (sample_number) {

        """
        mkdir -p reference
        cnvkit.py reference ${normal_coverage} \
            --fasta ${reference_fasta} \
            -o reference/normal_reference.cnn
        """

    } else {

        """
        mkdir -p reference
        cnvkit.py reference \
            --fasta ${reference_fasta} \
            -o reference/normal_reference.cnn \
            -t ${target} \
            -a ${antitarget}
        """
    }

    stub:
    """
    mkdir -p reference
    touch reference/normal_reference.cnn
    """

}