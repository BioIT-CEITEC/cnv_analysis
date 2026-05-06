process CLASSIFY_AND_ANNOTATE {

    tag "${meta.sample_name}"
    publishDir "structural_varcalls/${meta.sample_name}", mode: 'copy'
    conda "${moduleDir}/env.yaml" 

    input:
    tuple val(meta), path(merged_tsv), path(merged_bed)
    path annotation_tsv

    output:
    tuple val(meta), path("classified_and_annotated_CNVs/Scoresheet.txt"), path("classified_and_annotated_CNVs/*"), emit: cnvkit_annotated_calls

    script:

    def genomeBuild = params.assembly == "GRCh37" ? "hg19" : "hg38"

    """
    mkdir -p classified_and_annotated_CNVs
    python3 ${moduleDir}/ClassifyCNV.py --infile ${merged_bed} --outdir classified_and_annotated_CNVs --GenomeBuild ${genomeBuild} --precise
    Rscript ${moduleDir}/cnvAnnotateCNVkit.R ${merged_tsv} classified_and_annotated_CNVs/Scoresheet.txt ${annotation_tsv} classified_and_annotated_CNVs/${meta.sample_name}_final_CNVs_annotated.tsv classified_and_annotated_CNVs/${meta.sample_name}_final_CNVs_annotated.xlsx
    """


    stub:
    """
    mkdir -p classified_and_annotated_CNVs
    touch classified_and_annotated_CNVs/Scoresheet.txt
    touch classified_and_annotated_CNVs/${meta.sample_name}_final_CNVs_annotated.tsv
    """
}