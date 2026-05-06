process CNV_CALL_CNMOPS {

    tag "${meta.sample_name}"
    publishDir "structural_varcalls/${meta.sample_name}", mode: 'copy'
    conda "${moduleDir}/env.yaml"

    input:
    tuple val(meta), path(bam), path(bai)
    path cohort_data

    output:
    tuple val(meta), path("cnMOPS/${meta.sample_name}_cnMOPS.tsv"), emit: cnmops_cnvcalls

    script:

        """
        mkdir -p cnMOPS
        Rscript ${projectDir}/bin/cnMOPS_wrapper.R ${bam} ${cohort_data} ${meta.sample_name} cnMOPS/${meta.sample_name}_cnMOPS.tsv
        """

    stub:
        """
        mkdir -p cnMOPS
        touch cnMOPS/${meta.sample_name}_cnMOPS.tsv
        """

}