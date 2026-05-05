process COVERAGE_CALC {

    tag "${meta.sample_name}"
    conda "${moduleDir}/env.yaml"

    input:
    tuple val(meta), path(bam), path(bam_bai)
    path organism_reference
    path reference_index

    output:
    tuple val(meta), path("${meta.sample_name}.region_coverage.tsv"), emit: region_coverage

    script:
    """
    bedtools sort \
        -faidx ${reference_index} \
        -i ${organism_reference} \
        > ${meta.sample_name}.targets.sorted.bed

    SORT_ORDER=\$(samtools view -H ${bam} | awk '\$1=="@HD" {
        for (i=1; i<=NF; i++) {
            if (\$i ~ /^SO:/) {
                sub("SO:", "", \$i)
                print \$i
            }
        }
    }')

    if [ "\${SORT_ORDER}" = "coordinate" ]; then
        ln -s ${bam} ${meta.sample_name}.sorted.bam

        if [ -f "${bam_bai}" ]; then
            ln -s ${bam_bai} ${meta.sample_name}.sorted.bam.bai
        else
            samtools index ${meta.sample_name}.sorted.bam
        fi
    else
        samtools sort \
            -m 2G \
            -o ${meta.sample_name}.sorted.bam \
            ${bam}

        samtools index ${meta.sample_name}.sorted.bam
    fi

    bedtools coverage \
        -sorted \
        -a ${meta.sample_name}.targets.sorted.bed \
        -b ${meta.sample_name}.sorted.bam \
        -g ${reference_index} \
        > ${meta.sample_name}.region_coverage.tsv
    """

    stub:
    """
    touch ${meta.sample_name}.region_coverage.tsv
    """
}