process COHORT_ANALYSIS_CONIFER {
    conda "${moduleDir}/../env.yaml"

    input:
    path rpkm_files
    path probes

    output:
    path "conifer_results.h5", emit: conifer_cohort_analysis

    script:
    def conifer = "${moduleDir}/../conifer_v0.2.2/conifer.py"

    """
    mkdir -p rpkm_dir
    mv *.txt rpkm_dir/

    python ${conifer} \
        analyze \
        --probes ${probes} \
        --rpkm_dir rpkm_dir/ \
        --output conifer_results.h5 \
        --svd ${params.conifer_svd_components} \
        --min_rpkm ${params.conifer_min_rpkm}

    python -c "
import tables, sys
try:
    h5 = tables.open_file('conifer_results.h5', 'r')
    _ = h5.root.samples.samples
    h5.close()
except Exception as e:
    print '[ERROR] CoNIFER analyze produced incomplete output (missing /samples node): ' + str(e)
    sys.exit(1)
"
    """

    stub:
    """
    touch conifer_results.h5
    """
}