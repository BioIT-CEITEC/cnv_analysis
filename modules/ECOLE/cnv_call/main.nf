process CNV_CALL_ECOLE {
    publishDir "structural_varcalls/ECOLE", mode: 'copy'
    conda "${moduleDir}/../env.yaml"

    input:
    path bams
    path target_bed

    output:
    path("ecole_calls_output/*"), emit: ecole_cnvcalls

    script:

    """
    mkdir -p ecole_calls_output

    python ${moduleDir}/../ECOLE-0.2/scripts/ECOLE_call.py --model ecole \
    --input ./processed_samples \
    --output ./ecole_calls_output \
    --cnv exonlevel \
    --batch_size 16 \
    --normalize ${moduleDir}/../ECOLE-0.2/ecole_stats.txt \
    --gpu 0
    """

    stub:
    """
    mkdir -p ecole_calls_output
    for filename in ${bams}; do
        f=\$(basename -- "\$filename")
        touch "ecole_calls_output/\${f}.txt"
    done
    """

}