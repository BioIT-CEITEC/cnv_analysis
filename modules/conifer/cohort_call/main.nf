process COHORT_CALL_CONIFER {
    conda "${moduleDir}/../env.yaml"

    input:
    path cohort_hd5

    output:
    path "conifer_calls.txt", emit: conifer_cohort_calls

    script:
    """
    python ${moduleDir}/../conifer_v0.2.2/conifer.py \
        call \
        --input ${cohort_hd5} \
        --output conifer_calls.txt \
        --threshold ${params.conifer_threshold}
    """

    stub:
    """
    touch conifer_calls.txt
    """
}