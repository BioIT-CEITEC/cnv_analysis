process COHORT_ANALYSIS_CONIFER {
    conda "${moduleDir}/../env.yaml"

    input:
    path rpkm_files
    path probes

    output:
    path "conifer_results.h5", emit: conifer_cohort_analysis

    script:
    """
    python3 ${moduleDir}/conifer_v0.2.2/conifer.py \
        analyze \
        --probes ${probes} \
        --rpkm_dir ./ \
        --output conifer_results.h5 \
        --svd ${params.conifer_svd_components} \
        --min_rpkm ${params.conifer_min_rpkm}
    """

    stub:
    """
    touch conifer_probes.tsv
    """
}