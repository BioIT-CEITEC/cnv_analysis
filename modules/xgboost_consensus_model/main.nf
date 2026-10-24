process XGBOOST_CONSENSUS_SCORE {

    conda "${moduleDir}/env.yaml"
    tag "${meta.sample_name}"
    cache false

    // `sample_input` is either a single per-sample directory of caller output
    // (the standalone modelling/xgboost_pipeline entry point) or the flat list
    // of caller files the main pipeline carries in ch_all_varcalls_for_merging.
    // Both are reduced to one directory below; load_sample_calls() rglobs it,
    // so a flat directory of files and nested caller sub-directories are
    // equivalent as far as the model is concerned.
    input:
    tuple val(meta), path(sample_input)
    path model_json
    path panel_bed

    output:
    tuple val(meta), path("merged_variants/${meta.sample_name}_merged_target_consensus.tsv"),
                     path("merged_variants/${meta.sample_name}_merged_target_consensus.bed"), emit: final_tsvs
    tuple val(meta), path("merged_variants/${meta.sample_name}_smoothed_variants.tsv"),
                     path("merged_variants/${meta.sample_name}_smoothed_variants.bed"),       emit: merged_tsv
    tuple val(meta), path("merged_variants/${meta.sample_name}_per_exon_matrix.tsv"),         emit: raw_matrix
    tuple val(meta), path("merged_variants/${meta.sample_name}_all_scored_candidates.tsv"),   emit: all_scored
    tuple val(meta), path("merged_variants/${meta.sample_name}_refined_regions.tsv"),         emit: refined

    script:
    def sample = meta.sample_name
    def cutoff = params.xgb_cutoff
    def cutoff_arg = cutoff != null ? "--cutoff ${cutoff}" : ""

    // The complete scored table, always written, never filtered.
    def all_scored = "merged_variants/${sample}_all_scored_candidates.tsv"

    // Final region set: candidates broken down to supported panel exons and
    // re-smoothed. This is what gets annotated, so genes/classification/score
    // describe the coordinates they are attached to.
    def refined = "merged_variants/${sample}_refined_regions.tsv"

    // score_samples_xgb.py names its filtered output with the cutoff formatted
    // to 4 decimals; derive the same name here rather than globbing for it, so
    // a stale file from a different cutoff can never be picked up by mistake.
    // Locale.ROOT is required: Python's f"{c:.4f}" always emits a '.' decimal
    // separator, while a bare String.format follows the JVM default locale and
    // would emit '0,7499' under e.g. cs_CZ, silently missing the file.
    def filtered_tsv = cutoff != null
        ? "merged_variants/${sample}_smoothed_variants_filtered_${String.format(java.util.Locale.ROOT, '%.4f', cutoff as double)}.tsv"
        : ""
    """
    mkdir -p merged_variants

    # Reduce whichever input shape was staged to a single directory holding
    # this sample's caller files. A lone directory is used as-is (keeps caller
    # sub-directories, which load_sample_calls() prefers when several files
    # match the same CALLER_FILE_PATTERNS glob); loose files are symlinked into
    # sv_dir/<sample>/.
    sample_dir=sv_dir/${sample}
    mkdir -p "\$sample_dir"
    for f in ${sample_input}; do
        if [ -d "\$f" ]; then
            sample_dir="\$f"
        elif [ -f "\$f" ]; then
            ln -sf "\$(readlink -f "\$f")" "\$sample_dir"/
        fi
    done

    # Run model -- all scored candidates go to _smoothed_variants (pre-filter equivalent).
    # No fallback on failure: a sample with no calls already gets a header-only
    # table from score_samples_xgb.py itself, so a non-zero exit is a real error
    # (missing xgboost, unreadable model, parser bug) and must fail the task
    # rather than masquerade as a sample with no CNVs.
    python ${moduleDir}/score_samples_xgb.py \\
        --sample_dir  "\$sample_dir" \\
        --sample_name ${sample} \\
        --model       ${model_json} \\
        --out         merged_variants/${sample}_smoothed_variants.tsv \\
        ${cutoff_arg}

    # Preserve the complete scored candidate set under an explicit name BEFORE
    # the reshape at the end of this script overwrites _smoothed_variants.tsv
    # with the (possibly cutoff-filtered) table. This is the artefact to
    # re-threshold or audit against -- it must not depend on the incidental
    # ordering of the _per_exon_matrix copy below.
    cp merged_variants/${sample}_smoothed_variants.tsv ${all_scored}

    # Detect whether caller files use 'chr' prefix in chromosome names.
    # norm_chr() in cnv_consensus_model.py always normalises to 'chrN' form,
    # so we need to know the callers' original convention to restore it.
    # -L so the symlinks staged above are followed.
    strip_chr=1
    for f in \$(find -L "\$sample_dir" -type f); do
        chr_val=\$(awk '!/^#/ && NF>0 && NR>1 {print \$1; exit}' "\$f" 2>/dev/null | head -1)
        case "\$chr_val" in chr*) strip_chr=0; break;; esac
    done

    # _merged_target_consensus and the reshaped _smoothed_variants come from the
    # cutoff-filtered table when a cutoff was set, else from the complete scored
    # table. The existence guard matters even with a cutoff set: when a sample
    # yields no candidates at all, score_samples_xgb.py writes only a header to
    # --out and returns before writing any filtered file.
    src="${all_scored}"
    if [ -n "${filtered_tsv}" ]; then
        if [ -f "${filtered_tsv}" ]; then
            src="${filtered_tsv}"
        elif [ "\$(wc -l < ${all_scored})" -gt 1 ]; then
            # A cutoff was requested and there ARE candidates, so the filtered
            # file must exist. Missing means the derived name disagrees with
            # what score_samples_xgb.py wrote -- fail loudly rather than emit an
            # unfiltered table that looks filtered.
            echo "ERROR: expected filtered table ${filtered_tsv} not found" >&2
            ls -l merged_variants/ >&2
            exit 1
        fi
    fi

    # Turn the scored candidates into the FINAL region set BEFORE anything
    # annotates them: intersect with the capture BED, keep only exons carrying
    # real (non-bridging) caller support, re-smooth across small gaps. Union
    # candidates built by _merge_candidates can be inflated by a single
    # multi-megabase segmental call, and annotating that span attributes genes
    # and a classification to places no caller called.
    python ${moduleDir}/refine_regions.py \\
        --scored    "\$src" \\
        --sample_dir "\$sample_dir" \\
        --panel_bed ${panel_bed} \\
        --gap       ${params.smooth_gap_bp} \\
        --out       ${refined}

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
    cp ${all_scored} merged_variants/${sample}_per_exon_matrix.tsv

    # Reorder and normalise _smoothed_variants.tsv to match MERGE_VARIANT_CALLS schema exactly:
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
        ${refined} \\
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
    touch merged_variants/${meta.sample_name}_all_scored_candidates.tsv
    touch merged_variants/${meta.sample_name}_refined_regions.tsv
    """

}
