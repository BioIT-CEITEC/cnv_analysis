process REFERENCE_CNVKIT {

    publishDir "structural_varcalls/reference/cnvkit", mode: 'copy'

    input:
    path coverage_inputs // list of cnn files depending on calling_type and the number of samples
    tuple path(target), path(antitarget) // target and antitarget bed file produced in prepare_regions_cnvkit process
    path organism_fasta // fasta file of the organism genome

    output:
    path("reference/normal_reference.cnn"), emit: cnvkit_reference

    script:

    if (params.normal_coverage_inputs) {
        """
        mkdir -p reference
        cnvkit.py reference ${coverage_inputs} \
            --fasta ${organism_fasta} \
            -o reference/normal_reference.cnn
        """
    } else {
        """
        mkrdir -p reference
        cnvkit.py reference \
            --fasta ${organism_fasta} \
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