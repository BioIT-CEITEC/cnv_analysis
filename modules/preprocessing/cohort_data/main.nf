process COHORT_PREPROCESS {

    publishDir "cohort_data/"

    input:
    path cohort_tar

    output:
    path("cohort_data/cohort_info_tab.tsv"), emit: cohort_data

    script:
    """
    tar -xzf ${cohort_tar}
    """

    stub:
    """
    mkdir -p cohort_data
    touch cohort_data/cohort_info_tab.tsv
    """
}