process CNV_CALL_EXOMEDEPTH {

    tag "${meta.sample_name}"
    publishDir "structural_varcalls/${meta.sample_name}", mode: 'copy'
    conda "${moduleDir}/env.yaml"

    input:
    tuple val(meta), path(bam), path(bai)
    path cohort_data
    path reference_fasta

    output:
    tuple val(meta),path("exomeDepth/${meta.sample_name}_ExomeDepth.tsv"), emit: exomedepth_cnvcalls

    script:

        """
        mkdir -p exomeDepth

        Rscript ${projectDir}/bin/ExomeDepth_wrapper.R \
            ${bam} \
            ${cohort_data} \
            ${meta.sample_name} \
            ${reference_fasta} \
            exomeDepth/${meta.sample_name}_ExomeDepth.tsv

        """

    stub:
        """
        mkdir -p exomeDepth
        touch exomeDepth/${meta.sample_name}_ExomeDepth.tsv
        """

}