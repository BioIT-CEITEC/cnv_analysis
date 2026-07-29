process PER_REGION_COVERAGE_CALC {
    tag "${meta.sample_name}"
    publishDir "coverage_tracks/${meta.sample_name}", mode: 'copy', pattern: "*.region_coverage.tsv"
    conda "${moduleDir}/env.yaml"

    input:
    tuple val(meta), path(bam), path(bam_bai)
    path organism_bed
    path bed_annotation

    output:
    tuple val(meta), path("${meta.sample_name}.per_nucleotide_coverage.tsv"), emit: nucleotide_coverage
    tuple val(meta), path("${meta.sample_name}.region_coverage.tsv"),         emit: region_coverage

    script:
    """
    bedtools coverage \
        -d \
        -a ${organism_bed} \
        -b ${bam} \
        > ${meta.sample_name}.per_nucleotide_coverage.tsv

    awk 'BEGIN{OFS="\\t"}
    {
        n_bed = NF - 2
        key = \$1
        for (i = 2; i <= n_bed; i++) key = key OFS \$i
        if (!(key in seen)) {
            seen[key] = 1
            order[++nkeys] = key
            sum_d[key]   = 0
            sum_sq[key]  = 0
            n_cov[key]   = 0
            rlen[key]    = 0
        }
        sum_d[key]  += \$NF
        sum_sq[key] += \$NF * \$NF
        if (\$NF > 0) n_cov[key]++
        if (\$(NF-1) > rlen[key]) rlen[key] = \$(NF-1)
    }
    END {
        for (i = 1; i <= nkeys; i++) {
            k     = order[i]
            len   = rlen[k]
            cov   = n_cov[k] + 0
            mean  = (len > 0) ? sum_d[k] / len : 0
            var   = (len > 0) ? sum_sq[k] / len - mean * mean : 0
            stdev = (var > 0) ? sqrt(var) : 0
            frac  = (len > 0) ? cov / len : 0
            printf "%s\\t%.6f\\t%.6f\\t%d\\t%d\\t%.6f\\n", k, mean, stdev, cov, len, frac
        }
    }' ${meta.sample_name}.per_nucleotide_coverage.tsv \
        > ${meta.sample_name}.region_coverage.raw.tsv

    awk 'BEGIN{OFS="\\t"}
         NR==FNR { ann[\$1 "\\t" \$2 "\\t" \$3] = \$4 "\\t" \$5 "\\t" \$6 "\\t" \$7; next }
         { key = \$1 "\\t" \$2 "\\t" \$3
           print \$0 "\\t" (key in ann ? ann[key] : "NA\\tNA\\tNA\\tNA") }
    ' ${bed_annotation} ${meta.sample_name}.region_coverage.raw.tsv \
        > ${meta.sample_name}.region_coverage.tsv
    """

    stub:
    """
    touch ${meta.sample_name}.per_nucleotide_coverage.tsv
    touch ${meta.sample_name}.region_coverage.tsv
    """
}