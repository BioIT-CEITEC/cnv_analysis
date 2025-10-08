process SNP_AF_CALC {

    tag "${meta.sample_name}"
    conda "${moduleDir}/env.yaml"

    input:
    tuple val(meta), path(bam), path(bai)
    path reference_fasta
    path reference_index
    path reference_snps

    output:
    tuple val(meta), path("${meta.donor}_N.snpAF.tsv"), path("${meta.donor}_T.snpAF.tsv", optional: true), emit: snpAF

    script:

    //def hasNormals = params.panel_of_normals ?: false
    //def tumor_call = hasNormals ? "alleleCounter -r ${reference_fasta} -l ${reference_snps} -b ${tumor_bam} -o ${meta.donor}_T.snpAF.tsv" : "touch ${meta.donor}_T.snpAF.tsv"

    """
    alleleCounter -r ${reference_fasta} -l ${reference_snps} -b ${bam} -o ${meta.sample_name}.snpAF.tsv
    """
    stub:
    """
    touch ${meta.sample_name}.snpAF.tsv
    """
}