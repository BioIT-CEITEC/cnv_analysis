process CNV_CALL_PANELCNMOPS {

    tag "${meta.sample_name}"
    publishDir "structural_varcalls/panelcnMOPS", mode: 'copy'
    conda "${moduleDir}/env.yaml"

    input:
    tuple val(meta), path(bam), path(bai)
    path cohort_data

    output:
    tuple val(meta), path("cnv_calls/panelcnMOPS_CNV_${meta.sample_name}.tsv"), emit: panelcnMOPS_cnvcalls

    script:

        """
        mkdir -p cnv_calls

        Rscript ${projectDir}/bin/panelcnMOPS_wrapper.R \
        ${bam} \
        ${cohort_data} \
        ${meta.sample_name} \
        cnv_calls/panelcnMOPS_CNV_${meta.sample_name}.tsv

        """

    stub:
        """
        mkdir -p cohort_data
        touch cohort_data/target.bed
        touch cohort_data/antitarget.bed
        """

}