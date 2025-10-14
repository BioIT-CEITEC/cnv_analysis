process REFERENCE_CNVKIT {

    conda "${moduleDir}/../env.yaml"

    input:
    path reference_fasta // channel to the reference fasta file
    path normal_coverage // channel to the tumor coverage files
    path sample_coverage // channel to the normal coverage files
    tuple path(target), path(antitarget) // target and antitarget bed file produced in prepare_regions_cnvkit process
    val sample_number // if only tumor samples are provided at least 4 are required, otherwise will use antitarget.bed and targed.bed

    output:
    path("reference/normal_reference.cnn"), emit: cnvkit_reference

    script:

    def hasNormals = params.normal_tumor ?: false
    def coverage_files = hasNormals ? normal_coverage : sample_coverage

    if (hasNormals || sample_number) {

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