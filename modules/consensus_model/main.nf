process CONSENSUS_PREPARE_TRAINING_DIR {

    tag "${meta.sample_name}"

    input:
    tuple val(meta), path(varcall_files)

    output:
    path("${meta.sample_name}"), emit: sample_dir

    script:
    """
    # train_consensus_model.py expects one directory per sample, named after the
    # sample, holding that sample's caller files (matched by CALLER_FILE_PATTERNS).
    mkdir -p ${meta.sample_name}
    for f in ${varcall_files}; do
        [ -f "\$f" ] && cp -L "\$f" ${meta.sample_name}/
    done
    """

    stub:
    """
    mkdir -p ${meta.sample_name}
    """

}

process CONSENSUS_MODEL_TRAIN {

    conda "${moduleDir}/env.yaml"
    publishDir "consensus_model", mode: 'copy'

    input:
    path sample_dirs, stageAs: 'sv_dir/*'
    path gt_dir

    output:
    path("cnv_consensus_model.pkl"),  emit: model
    path("cnv_consensus_model.json"), emit: model_json
    path("cnv_train_*.tsv"),          emit: training_stats

    script:
    def callers_arg = params.consensus_train_callers ? "--callers ${params.consensus_train_callers}" : ""
    def cutoff_arg  = params.consensus_cutoff != null ? "--score_cutoff ${params.consensus_cutoff}" : ""
    """
    python ${moduleDir}/train_consensus_model.py \\
        --sim_dir     sv_dir \\
        --gt_dir      ${gt_dir} \\
        --out_model   cnv_consensus_model.json \\
        --out_prefix  cnv_train \\
        --weight_cap  ${params.consensus_train_weight_cap} \\
        --cap_mode    ${params.consensus_train_cap_mode} \\
        --l2          ${params.consensus_train_l2} \\
        --test_size   ${params.consensus_train_test_size} \\
        --random_seed ${params.consensus_train_seed} \\
        ${callers_arg} \\
        ${cutoff_arg}
    """

    stub:
    """
    touch cnv_consensus_model.pkl
    touch cnv_consensus_model.json
    touch cnv_train_caller_stats.tsv
    """

}

process CONSENSUS_MODEL_SCORE {

    conda "${moduleDir}/env.yaml"
    tag "${meta.sample_name}"
    cache false

    input:
    tuple val(meta), path(varcall_files)
    path model_pkl

    output:
    tuple val(meta), path("merged_variants/${meta.sample_name}_merged_target_consensus.tsv"),
                     path("merged_variants/${meta.sample_name}_merged_target_consensus.bed"), emit: final_tsvs
    tuple val(meta), path("merged_variants/${meta.sample_name}_smoothed_variants.tsv"),
                     path("merged_variants/${meta.sample_name}_smoothed_variants.bed"),       emit: merged_tsv
    tuple val(meta), path("merged_variants/${meta.sample_name}_per_exon_matrix.tsv"),         emit: raw_matrix

    script:
    def sample    = meta.sample_name
    def cutoff_arg = params.consensus_cutoff != null ? "--cutoff ${params.consensus_cutoff}" : ""
    """
    # Stage caller files into the per-sample directory score_samples.py expects
    mkdir -p sv_dir/${sample}
    for f in ${varcall_files}; do
        [ -f "\$f" ] && ln -sf \$(readlink -f "\$f") sv_dir/${sample}/
    done

    mkdir -p merged_variants

    # Run model — all scored candidates go to _smoothed_variants (pre-filter equivalent)
    python ${moduleDir}/score_samples.py \\
        --sv_dir sv_dir \\
        --model  ${model_pkl} \\
        --out    merged_variants/${sample}_smoothed_variants.tsv \\
        ${cutoff_arg} \\
    || {
        # No calls found: write headers-only file so outputs always exist
        printf 'sample\\tchr\\tstart\\tend\\ttype\\tn_callers\\tcallers_supporting\\tconsensus_score\\tcn_label\\n' \\
            > merged_variants/${sample}_smoothed_variants.tsv
    }

    # Detect whether caller files use 'chr' prefix in chromosome names.
    # norm_chr() in cnv_consensus_model.py always normalises to 'chrN' form,
    # so we need to know the callers' original convention to restore it.
    strip_chr=1
    for f in ${varcall_files}; do
        [ -f "\$f" ] || continue
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

    # Full TSV — all columns preserved, only chr prefix stripped if callers don't use it.
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

    # _per_exon_matrix: no direct equivalent — reuse all-scored file for downstream compatibility
    cp merged_variants/${sample}_smoothed_variants.tsv \\
       merged_variants/${sample}_per_exon_matrix.tsv

    # Reorder and normalise _smoothed_variants.tsv to match MERGE_VARIANT_CALLS schema exactly:
    #   sample(1) CHROM(2) consensus_type(3) START(4) END(5) n_targets(6) n_callers(7)
    #   callers(8) target_names(9) BED_gene_name(10) gene_biotype(11) cn_label(12) consensus_score(13)
    #
    # Why sample first: cnvAnnotateCNVkit.R detects UCSC/ENS by inputDF[1,1]; when sample is
    # first the name ("BR-2604") does not start with "chr" → UCSC mode. UCSC mode strips "chr"
    # from ClassifyCNV VariantIDs and filters the GTF TSV with Ensembl chromosome names ("1","2"…),
    # which matches the GRCh38.tsv used in the pipeline. ENS mode would filter by "chr1","chr2"…
    # and find nothing in an Ensembl-format GTF TSV → foverlaps produces 0 rows → empty output.
    #
    # chr stripping: norm_chr() in cnv_consensus_model.py always emits "chrN"; strip it here so
    # chromosome names match the Ensembl-format BED/GTF used everywhere else in the pipeline.
    # callers_supporting uses "," separator; convert to ";" to match merging_and_smoothing.R.
    # n_targets, target_names, BED_gene_name, gene_biotype filled NA (not from model).
    # cn_label populated from caller CN columns (panelcnMOPS CN integer, XHMM CN, gatk CN; NA for exomeDepth/conifer).
    # consensus_score kept as extra column at position 13 (ignored by combine_final_tables.R).
    #
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
    # Use 1-based START so VariantID from ClassifyCNV (chr_start_end_type) and from
    # paste(CHR,START,STOP,CNVtype) in cnvAnnotateCNVkit.R are identical after UCSC chr-stripping.
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
