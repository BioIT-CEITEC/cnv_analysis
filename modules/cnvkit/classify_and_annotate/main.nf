process CNVKIT_CLASSIFY_AND_ANNOTATE {

    tag "${meta.donor}"
    publishDir "structural_varcalls/${meta.donor}/cnvkit", mode: 'copy'

    conda "${moduleDir}/../env.yaml" 

    input:
    tuple val(meta), path(cns_calls), path(preannot_bed), path(bed_calls), path(vcf_calls), path(preannot_tsv)
    path annotation_tsv

    output:
    tuple val(meta), path("classify/*"), path("annotate/*"), emit: cnvkit_annotated_calls

    script:

    def genomeBuild = params.assembly == "GRCh37" ? "hg19" : "hg38"

    """
    mkdir -p annotate
    python3 ${moduleDir}/ClassifyCNV.py --infile ${preannot_bed} --outdir classify --GenomeBuild ${genomeBuild} --precise
    Rscript ${moduleDir}/cnvAnnotateCNVkit.R ${preannot_tsv} classify/Scoresheet.txt ${annotation_tsv} annotate/${meta.donor}_final_CNVs_annotated.tsv annotate/${meta.donor}_final_CNVs_annotated.xlsx
    """

    stub:
    """
    mkdir -p call
    touch call/CNV_calls.cns
    """
}