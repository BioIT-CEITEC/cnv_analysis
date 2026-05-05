process EXTRACT_SAMPLE_FREEC {
    tag "${meta.sample_name}"
    publishDir "structural_varcalls/${meta.sample_name}/freec", mode: 'copy'
    conda "${moduleDir}/../env.yaml"

    input:
    tuple val(meta), path(cnvs), path(config), path(ratio)

    output:
    tuple val(meta), path("${meta.sample_name}_freec.tsv"), emit: cnv_tsv

    script:
    """
    python3 << CODE
rows = []
with open("${cnvs}") as fh:
    for line in fh:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split("\\t")
        if len(fields) < 5:
            continue
        chrom, start, end, cn, status = fields[0], fields[1], fields[2], fields[3], fields[4].strip().lower()
        if status == "gain":
            cnv_type = "duplication"
        elif status == "loss":
            cnv_type = "deletion"
        else:
            continue
        rows.append((chrom, start, end, cn, cnv_type))

with open("${meta.sample_name}_freec.tsv", "w") as fh:
    fh.write("chromosome\\tstart\\tend\\tCN\\ttype\\n")
    for r in rows:
        fh.write("\\t".join(r) + "\\n")
CODE
    """

    stub:
    """
    touch ${meta.sample_name}_freec.tsv
    """
}
