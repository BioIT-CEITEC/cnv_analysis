process GET_COVERAGE_CNVKIT {
    tag "${meta.sample_name}"

    conda "${moduleDir}/../env.yaml"

    input:
    tuple val(meta), path(bam), path(bam_bai)
    tuple path(target), path(antitarget)

    output:
    tuple val(meta), 
          path("coverage/${meta.sample_name}.targetcoverage.cnn"),
          path("coverage/${meta.sample_name}.antitargetcoverage.cnn"),
          emit: coverage 

    script:

    """
    mkdir -p coverage
    cnvkit.py coverage ${bam} ${target} -o coverage/${meta.sample_name}.targetcoverage.cnn
    cnvkit.py coverage ${bam} ${antitarget} -o coverage/${meta.sample_name}.antitargetcoverage.cnn
    """

    stub:
    """
    mkdir -p coverage
    touch coverage/${meta.sample_name}.targetcoverage.cnn
    touch coverage/${meta.sample_name}.antitargetcoverage.cnn
    """
}