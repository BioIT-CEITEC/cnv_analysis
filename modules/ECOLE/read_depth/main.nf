process READ_DEPTH_ECOLE {

    conda "${moduleDir}/../env.yaml"

    input:
    path bams
    path target_bed

    output:
    path("processed_samples/"), emit: ecole_preprocessed_samples

    script:

    """
    mkdir -p read_depths
    mkdir -p processed_samples

    for filename in ./*.bam; do
        f=\$(basename -- "\$filename")
        sambamba depth base -L ${target_bed} "\$filename" > "read_depths/\${f}.txt"
    done

    python ${moduleDir}/../ECOLE-0.2/scripts/preprocess_sample.py \
        --readdepth ./read_depths \
        --output ./processed_samples --target ${target_bed}
    """

    stub:
    """
    mkdir -p read_depths
    for filename in ./*.bam; do
        f=\$(basename -- "\$filename")
        touch "read_depths/\${f}.txt"
    done
    """

}