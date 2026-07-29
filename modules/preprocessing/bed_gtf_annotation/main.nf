process BED_GTF_ANNOTATION {

    publishDir "structural_varcalls/bed_annotation", mode: 'copy'
    conda "${moduleDir}/env.yaml"

    input:
    path bed_file
    path gtf_file

    output:
    path "bed_gtf_annotation.tsv", emit: annotation

    script:
    def annotate = params.organism_gtf ? true : false
    def gtf_cmd  = annotate
        ? "python3 ${projectDir}/bin/build_bed_gtf_annotation.py ${gtf_file} ${bed_file} bed_gtf_annotation.tsv"
        : "awk 'BEGIN{OFS=\"\\t\"} !/^#/ && NF>0 {print \$1,\$2,\$3,\"NA\",\"NA\",\"NA\",\"NA\"}' ${bed_file} > bed_gtf_annotation.tsv"
    """
    ${gtf_cmd}
    """

    stub:
    """
    awk 'BEGIN{OFS="\\t"} !/^#/ && NF>0 {print \$1,\$2,\$3,"NA","NA","NA","NA"}' ${bed_file} \
        > bed_gtf_annotation.tsv
    """
}
