process CNV_CALL_EXOMEDEPTH {

    tag "${meta.sample_name}"
    publishDir "structural_varcalls/exomeDepth", mode: 'copy'
    conda "${moduleDir}/env.yaml"

    input:
    tuple val(meta), path(bam), path(bai)
    path cohort_data
    path reference_fasta

    output:
    path("cnv_calls/ExomeDepth_CNV_${meta.sample_name}.tsv"), emit: exomedepth_cnvcalls

    script:

        """
        mkdir -p cnv_calls

        Rscript ${projectDir}/bin/ExomeDepth_wrapper.R \
            ${bam} \
            ${cohort_data} \
            ${meta.sample_name} \
            ${reference_fasta} \
            cnv_calls/ExomeDepth_CNV_${meta.sample_name}.tsv

        """

    stub:
        """
        mkdir -p cohort_data
        touch cohort_data/target.bed
        touch cohort_data/antitarget.bed
        """

}