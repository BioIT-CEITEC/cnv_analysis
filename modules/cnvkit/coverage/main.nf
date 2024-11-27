process GET_COVERAGE_CNVKIT {

    tag "${meta.id}"
    publishDir "structural_varcalls/${meta.id}/cnvkit", mode: 'copy'

    input:
    tuple val(meta), path(bam_tumor), path(bam_bai_tumor), path(bam_normal), path(bam_bai_normal) // mandatory [ [meta],[bam_tumor],[bam_bai_tumor],[bam_normal],[bam_bai_normal] ]
    tuple path(target), path(antitarget) // target and antitarget bed file produced in prepare_regions_cnvkit process

    output:
    tuple val(meta), 
    path("coverage/tumor.targetcoverage.cnn"),
    path("coverage/tumor.antitargetcoverage.cnn"),
    path("coverage/normal.targetcoverage.cnn", optional: true),
    path("coverage/normal.antitargetcoverage.cnn", optional: true),
    emit: coverage 

    // this optional comes handy here but it has to be handled in the the subworkflow as we need the to set the last two as empty later

    script:

    def normal_cmd = params.normal_tumor ? "cnvkit.py coverage ${bam_normal} ${target} -o coverage/normal.targetcoverage.cnn && cnvkit.py coverage ${bam_normal} ${antitarget} -o coverage/normal.antitargetcoverage.cnn" : ""

    """
    mkdir -p coverage
    cnvkit.py coverage ${bam_tumor} ${target} -o coverage/tumor.targetcoverage.cnn
    cnvkit.py coverage ${bam_tumor} ${antitarget} -o coverage/tumor.antitargetcoverage.cnn
    ${normal_cmd}
    """

    stub:
    """
    mkdir -p coverage
    touch coverage/tumor.targetcoverage.cnn
    touch coverage/tumor.antitargetcoverage.cnn
    touch coverage/normal.targetcoverage.cnn
    touch coverage/normal.antitargetcoverage.cnn
    """
}