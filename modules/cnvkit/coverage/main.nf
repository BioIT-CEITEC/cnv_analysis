process GET_COVERAGE_CNVKIT {
    tag "${meta.sample_name}"
    publishDir "structural_varcalls/$meta.sample_name/cnvkit", mode: 'copy'

    conda "${moduleDir}/../env.yaml"

    input:
    tuple val(meta), path(bam), path(bam_bai)
    tuple path(target), path(antitarget)

    output:
    tuple val(meta), 
          path("coverage/${meta.sample_name}.tumor.targetcoverage.cnn"), 
          path("coverage/${meta.sample_name}.tumor.antitargetcoverage.cnn"), 
          path("coverage/${meta.sample_name}.normal.targetcoverage.cnn"),
          path("coverage/${meta.sample_name}.normal.antitargetcoverage.cnn"),
          emit: coverage 

    script:

    """
    mkdir -p coverage
    cnvkit.py coverage ${bam} ${target} -o coverage/${meta.sample_name}.tumor.targetcoverage.cnn
    cnvkit.py coverage ${bam} ${antitarget} -o coverage/${meta.sample_name}.tumor.antitargetcoverage.cnn
    """

    stub:
    """
    mkdir -p coverage
    touch coverage/${meta.sample_name}.tumor.targetcoverage.cnn
    touch coverage/${meta.sample_name}.tumor.antitargetcoverage.cnn
    touch coverage/${meta.sample_name}.normal.targetcoverage.cnn
    touch coverage/${meta.sample_name}.normal.antitargetcoverage.cnn
    """
}