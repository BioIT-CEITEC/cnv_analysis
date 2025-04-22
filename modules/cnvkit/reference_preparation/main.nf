process REFERENCE_CNVKIT {

    publishDir "structural_varcalls/reference/cnvkit", mode: 'copy'

    input:
    path reference_fasta // channel to the reference fasta file
    path tumor_coverage // channel to the tumor coverage files
    path normal_coverage // channel to the normal coverage files
    tuple path(target), path(antitarget) // target and antitarget bed file produced in prepare_regions_cnvkit process
    val sample_number // if only tumor samples are provided at least 4 are required, otherwise will use antitarget.bed and targed.bed

    output:
    path("reference/normal_reference.cnn"), emit: cnvkit_reference

    script:

    def coverage_files = params.tumor_normal ? normal_coverage : tumor_coverage

    if (params.tumor_normal || sample_number) {

        """
        mkdir -p reference
        cnvkit.py reference $coverage_files \
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