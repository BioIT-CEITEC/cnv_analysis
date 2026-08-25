process XGBOOST_CONSENSUS_SCORE {

    conda "${moduleDir}/env.yaml"
    tag "${meta.sample_name}"
    cache false

    input:
    tuple val(meta), path(sample_dir)
    path model_json
    // Staged into the task's working directory so score_samples_xgb.py can
    // import it via `sys.path.insert(0, ".")` -- the script never hardcodes
    // where cnv_consensus_model.py physically lives on disk.
    path cnv_consensus_model_py

    output:
    tuple val(meta), path("merged_variants/${meta.sample_name}_merged_target_consensus.tsv"),
                     path("merged_variants/${meta.sample_name}_merged_target_consensus.bed"), emit: final_tsvs
    tuple val(meta), path("merged_variants/${meta.sample_name}_smoothed_variants.tsv"),
                     path("merged_variants/${meta.sample_name}_smoothed_variants.bed"),       emit: merged_tsv
    tuple val(meta), path("merged_variants/${meta.sample_name}_per_exon_matrix.tsv"),         emit: raw_matrix

    script:
    def sample = meta.sample_name
    def cutoff_arg = params.xgb_cutoff != null ? "--cutoff ${params.xgb_cutoff}" : ""
    """
    mkdir -p merged_variants

    # Run model -- all scored candidates go to _smoothed_variants (pre-filter equivalent)
    python ${moduleDir}/score_samples_xgb.py \\
        --sample_dir  ${sample_dir} \\
        --sample_name ${sample} \\
        --model       ${model_json} \\
        --out         merged_variants/${sample}_smoothed_variants.tsv \\
        ${cutoff_arg} \\
    || {
        # No calls found: write headers-only file so outputs always exist
        printf 'sample\\tchr\\tstart\\tend\\ttype\\tn_callers\\tcallers_supporting\\tconsensus_score\\tcn_label\\n' \\
            > merged_variants/${sample}_smoothed_variants.tsv
    }

    # Detect whether caller files use 'chr' prefix in chromosome names.
    # norm_chr() in cnv_consensus_model.py always normalises to 'chrN' form,
    # so we need to know the callers' original convention to restore it.
    # (Same probe as CONSENSUS_MODEL_SCORE in modules/consensus_model/main.nf.)
    strip_chr=1
    for f in \$(find ${sample_dir} -maxdepth 2 -type f); do
        chr_val=\$(awk '!/^#/ && NF>0 && NR>1 {print \$1; exit}' "\$f" 2>/dev/null | head -1)
        case "\$chr_val" in chr*) strip_chr=0; break;; esac
    done

    # _merged_target_consensus: filtered (above cutoff) when a cutoff was set, else all calls
    filtered=\$(ls merged_variants/${sample}_smoothed_variants_filtered_*.tsv 2>/dev/null | head -1)
    if [ -n "\$filtered" ]; then
        src="\$filtered"
    else
        src="merged_variants/${sample}_smoothed_variants.tsv"
    fi

    # Full TSV -- all columns preserved, only chr prefix stripped if callers don't use it.
    # TSV source columns: sample(1) chr(2) start(3) end(4) type(5) n_callers(6) callers_supporting(7) consensus_score(8)
    awk -v strip="\$strip_chr" 'BEGIN {OFS="\\t"}
        NR == 1 { print; next }
        {
            chrom = \$2
            if (strip && substr(chrom,1,3) == "chr") chrom = substr(chrom,4)
            \$2 = chrom
            print
        }' "\$src" > merged_variants/${sample}_merged_target_consensus.tsv

    # Matching BED (0-based start)
    awk -v strip="\$strip_chr" 'BEGIN {OFS="\\t"}
        NR > 1 {
            chrom = \$2
            if (strip && substr(chrom,1,3) == "chr") chrom = substr(chrom,4)
            print chrom, \$3 - 1, \$4, \$5
        }' "\$src" > merged_variants/${sample}_merged_target_consensus.bed

    # _per_exon_matrix: no direct equivalent -- reuse all-scored file for downstream compatibility
    cp merged_variants/${sample}_smoothed_variants.tsv \\
       merged_variants/${sample}_per_exon_matrix.tsv

    # Reorder and normalise _smoothed_variants.tsv to match MERGE_VARIANT_CALLS schema exactly
    # (identical reshape to CONSENSUS_MODEL_SCORE, modules/consensus_model/main.nf):
    #   sample(1) CHROM(2) consensus_type(3) START(4) END(5) n_targets(6) n_callers(7)
    #   callers(8) target_names(9) BED_gene_name(10) gene_biotype(11) cn_label(12) consensus_score(13)
    # Source col order: sample(1) chr(2) start(3) end(4) type(5) n_callers(6) callers_supporting(7) consensus_score(8) cn_label(9)
    awk 'BEGIN {OFS="\\t"}
        NR == 1 { print "sample","CHROM","consensus_type","START","END","n_targets","n_callers","callers","target_names","BED_gene_name","gene_biotype","cn_label","consensus_score"; next }
        {
            chrom = \$2; sub(/^chr/, "", chrom)
            callers = \$7; gsub(",", ";", callers)
            cn_label = (NF >= 9) ? \$9 : "NA"
            print \$1, chrom, \$5, \$3, \$4, "NA", \$6, callers, "NA", "NA", "NA", cn_label, \$8
        }' \\
        "\$src" \\
        > merged_variants/${sample}_smoothed_variants_renamed.tsv
    mv merged_variants/${sample}_smoothed_variants_renamed.tsv \\
       merged_variants/${sample}_smoothed_variants.tsv

    # BED for _smoothed_variants: after reorder col 2=CHROM, col 3=consensus_type, col 4=START, col 5=END
    awk 'BEGIN {OFS="\\t"} NR > 1 {print \$2, \$4, \$5, \$3}' \\
        merged_variants/${sample}_smoothed_variants.tsv \\
        > merged_variants/${sample}_smoothed_variants.bed
    """

    stub:
    """
    mkdir -p merged_variants
    touch merged_variants/${meta.sample_name}_merged_target_consensus.tsv
    touch merged_variants/${meta.sample_name}_merged_target_consensus.bed
    touch merged_variants/${meta.sample_name}_smoothed_variants.tsv
    touch merged_variants/${meta.sample_name}_smoothed_variants.bed
    touch merged_variants/${meta.sample_name}_per_exon_matrix.tsv
    """

}
