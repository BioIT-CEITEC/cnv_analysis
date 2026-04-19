process CLASSIFY_AND_ANNOTATE {

    tag "${meta.sample_name}"
    publishDir "annotated_varcalls/", mode: 'copy'
    conda "${moduleDir}/env.yaml" 

    input:
    tuple val(meta), path(merged_tsv), path(merged_bed)
    path annotation_tsv

    output:
    tuple val(meta), path("classify/*"), path("annotate/*"), emit: cnvkit_annotated_calls

    script:

    def genomeBuild = params.assembly == "GRCh37" ? "hg19" : "hg38"

    """
    mkdir -p annotate
    python3 ${moduleDir}/ClassifyCNV.py --infile ${merged_bed} --outdir classify --GenomeBuild ${genomeBuild} --precise
    Rscript ${moduleDir}/cnvAnnotateCNVkit.R ${merged_tsv} classify/Scoresheet.txt ${annotation_tsv} annotate/${meta.sample_name}_final_CNVs_annotated.tsv annotate/${meta.sample_name}_final_CNVs_annotated.xlsx
    """


    stub:
    """
    mkdir -p call
    touch call/CNV_calls.cns
    """
}