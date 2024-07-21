process COHORT_PREPROCESS {

    publishDir "cohort_data/cohort_data/jabConTool"

    input:
    path(cohort_tar)

    output:
    path("cohort_info_tab.tsv"), emit: cohort_data

    script:

    """
    tar -xzf $cohort_tar
    """
}