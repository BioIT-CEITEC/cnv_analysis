process CNV_CALL_ECOLE {
    publishDir "structural_varcalls/all_samples", mode: 'copy'
    conda "${moduleDir}/../env.yaml"

    input:
    path bams
    path target_bed

    output:
    path("ECOLE/*"), emit: ecole_varcalls

    script:
    def ecole_path = "${moduleDir}/../ECOLE-0.2/scripts/ECOLE_call.py"
    """
    mkdir -p ECOLE
    python ${ecole_path} --model ecole \
    --input ./processed_samples \
    --output ./ECOLE \
    --cnv exonlevel \
    --batch_size 16 \
    --normalize ${moduleDir}/../ECOLE-0.2/ecole_stats.txt \
    --gpu 0
    """

    stub:
    """
    mkdir -p ECOLE
    for filename in ${bams}; do
        f=\$(basename -- "\$filename")
        touch "ECOLE/\${f}.txt"
    done
    """

}