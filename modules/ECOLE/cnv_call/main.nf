process CNV_CALL_ECOLE {
    publishDir "structural_varcalls/ECOLE", mode: 'copy'
    conda "${moduleDir}/env.yaml"

    input:
    path bams
    path target_bed

    output:
    path("ecole_calls_output/*"), emit: ecole_cnvcalls

    script:

    """
    mkdir -p read_depths

    python ${moduleDir}/ECOLE-02/scripts/ECOLE_call.py --model ecole --input ./processed_samples --output ./ecole_calls_output --cnv merged --batch_size 16 --normalize ecole_stats.txt --gpu 0
    """

    stub:
    """
    mkdir -p read_depths
    for filename in ${bams}; do
        f=\$(basename -- "\$filename")
        touch "read_depths/\${f}.txt"
    done
    """

}