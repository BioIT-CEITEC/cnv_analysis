process CNV_CALL_CNMOPS {

    tag "${meta.sample_name}"
    publishDir "structural_varcalls/cnMOPS", mode: 'copy'
    conda "${moduleDir}/env.yaml"

    input:
    tuple val(meta), path(bam), path(bai)
    path cohort_data

    output:
    path("cnv_calls/cnMOPS_CNV_${meta.sample_name}.tsv"), emit: cnmops_cnvcalls

    script:

        """
        mkdir -p cnv_calls

        Rscript ${projectDir}/bin/cnMOPS_wrapper.R ${bam} ${cohort_data} ${meta.sample_name} cnv_calls/cnMOPS_CNV_${meta.sample_name}.tsv

        """

    stub:
        """
        mkdir -p cohort_data
        touch cohort_data/target.bed
        touch cohort_data/antitarget.bed
        """

}