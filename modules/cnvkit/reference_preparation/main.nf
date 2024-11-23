process REFERENCE_CNVKIT {

    publishDir "structural_varcalls/reference/cnvkit", mode: 'copy'

    input:
    tuple val(meta), path(bam_tumor), path(bam_bai_tumor), path(bam_normal), path(bam_bai_normal)
    tuple path(target), path(antitarget)
    path reference

    output:
    path("normal_reference.cnn"), emit: reference

    script:

    def normal_ref = params.normal_tumor ? "${bam_normal}" : ""

    if (params.normal_coverage_inputs) {
        """
        cnvkit.py reference ${bam_tumor} ${normal_ref} \
            --fasta ${reference} \
            -o normal_reference.cnn
        """
    } else {
        """
        cnvkit.py reference \
            --fasta ${reference} \
            -o normal_reference.cnn \
            -t ${target} \
            -a ${antitarget}
        """
    }

    stub:
    """
    touch normal_reference.cnn
    """

}