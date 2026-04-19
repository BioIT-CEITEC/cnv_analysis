process MISMATCH_PILEUP {
    tag "${meta.sample_name}"
    publishDir "pseudogene/mismatch_pileup", mode: 'copy'
    conda "${moduleDir}/../env.yaml"

    input:
    tuple val(meta), path(realigned_dir), path(bam), path(bai)
    path diff_dir
    path region_bed

    output:
    tuple val(meta), path("${meta.sample_name}_mismatch_pileup.tsv"), emit: pileup_tsv

    script:
    """
    python3 ${projectDir}/bin/mismatch_pileup.py \
        --diff_dir    ${diff_dir} \
        --realigned_dir ${realigned_dir} \
        --bed         ${region_bed} \
        --sample      ${meta.sample_name} \
        --output      ${meta.sample_name}_mismatch_pileup.tsv \
        
    """
}
