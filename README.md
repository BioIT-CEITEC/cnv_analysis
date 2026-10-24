# cnv_analysis — germline CNV/SV calling pipeline

A Nextflow (DSL2) pipeline that runs a panel of independent CNV callers over targeted
(panel/WES) BAM files, combines their calls into a per-sample consensus, annotates and
classifies the result, and emits cohort-level tables plus an HTML report.

Two combination strategies are available:

| Strategy | Enabled by | Implementation |
|---|---|---|
| Rule-based merge & smooth | `use_consensus = false` | `bin/merging_and_smoothing.R` via `MERGE_VARIANT_CALLS` |
| **XGBoost consensus model** (default in `nextflow.config`) | `use_consensus = true` | `modules/xgboost_consensus_model/` — tuned gradient-boosted trees over 6 callers, trained offline (see [Consensus model](#consensus-model)) |

Both write the same output files, so annotation and the cohort tables are unaffected by the
choice.

---

## Contents

- [Requirements](#requirements)
- [Repository layout](#repository-layout)
- [Quick start](#quick-start)
- [Input data](#input-data)
- [Configuration](#configuration)
- [What the pipeline does](#what-the-pipeline-does)
- [Output layout](#output-layout)
- [Consensus model](#consensus-model)
  - [Candidate regions](#candidate-regions)
  - [Features and model](#features-and-model)
  - [Scoring (`use_consensus = true`)](#scoring-use_consensus--true)
  - [Training the model](#training-the-model)
  - [Pitfalls and caveats](#pitfalls-and-caveats)
- [Cohort data reuse](#cohort-data-reuse)
- [Resources and executor tuning](#resources-and-executor-tuning)

---

## Requirements

- **Nextflow** ≥ 22.10 — DSL2
- **Conda** (or Mamba) — every process declares its own `env.yaml`; Nextflow creates the
  environments on first run (`conda.createTimeout = '2 h'` is set because some are large)
- A CUDA GPU is **optional**; only ECOLE benefits from one (`--gpu 0`), everything else is CPU
- Disk: caller environments plus per-sample read-depth intermediates. ECOLE's environment is a
  full Anaconda export pinned to Python 3.8 / CUDA 10.1 and is by far the heaviest.

## Repository layout

```
main.nf                     entry point → workflow CNV_ANALYSIS → PANEL_WES
nextflow.config             all params, conda/executor settings, per-process resources
local.config                site overrides (tmp dirs, conda cacheDir, scratch)
workflow.config.json        descriptor for the institutional workflow manager (LOG)
workflows/
  PANEL.nf                  the panel/WES workflow — reference loading, caller toggles, wiring
subworkflows/<tool>/main.nf per-caller orchestration
modules/<tool>/main.nf      individual processes (+ env.yaml, helper scripts)
  xgboost_consensus_model/  the consensus model — scorer, shipped model, and its
                            offline training/evaluation code under modelling/
bin/                        R/Python helper scripts called by the modules
lib/Utils.groovy            sample-sheet parsing and the varcall channel helpers
```

## Quick start

```bash
# 1. put indexed BAMs in <repo>/mapped/  (see "Input data")
# 2. edit nextflow.config: references, lib_ROI, samples, caller toggles
# 3. run
nextflow run main.nf -c local.config -resume
```

Useful variants:

```bash
nextflow run main.nf -stub -resume        # exercise the channel topology; processes that define a
                                          # stub: block (73 of 76) become no-ops

nextflow run main.nf --use_ecole false    # override any param on the command line
nextflow run main.nf -with-report -with-trace
```

## Input data

Samples are declared in the `samples` map in `nextflow.config`; `Utils.parseInputVC()`
resolves each entry to a fixed location:

```
${projectDir}/mapped/${sample_name}.bam
${projectDir}/mapped/${sample_name}.bam.bai
```

`sample_name` must therefore be the **full BAM basename without `.bam`**, and the BAMs (or
symlinks to them) must live under `mapped/` inside the pipeline directory:

```groovy
samples = [
    '007_CZE3xBRCA1328_run137': [sample_name: '007_CZE3xBRCA1328_run137.picard.sorted.RG.rmdup'],
    ...
]
```

The map key is only an index label; `sample_name` is what propagates into `meta.sample_name`
and hence into every output filename. Remember this when preparing ground truth for the
consensus model — see [Pitfalls and caveats](#pitfalls-and-caveats).

## Configuration

All parameters live in the `params { }` block of `nextflow.config`.

**References**

| Param | Purpose |
|---|---|
| `assembly` | e.g. `hg19` — informational/annotation selection |
| `organism_fasta`, `organism_dict`, `organism_fasta.fai` (derived) | genome sequence and dictionary |
| `organism_ploidy_priors` | expected copies per chromosome (GATK gCNV) |
| `organism_gtf`, `organism_gtf_tsv` | annotation; some tools need the tabular form |
| `organism_cytoband` | cytoband table (annotation) |
| `organism_dna_panel` | the panel BED (regions of interest) |
| `organism_snps_panel` | known SNPs inside the ROI (jabCoNtool allele fractions) |
| `lib_ROI` | name of the panel region list, must match the BED name |

**Caller toggles** — each is an independent `true`/`false`:

`use_cnvkit`, `use_jabcontool`, `use_gatk`, `use_ecole`, `use_xhmm`, `use_conifer`,
`use_controlfreec`, `use_panelcnmops`, `use_cnmops`, `use_exomedepth`

**Combination**: `use_consensus`, plus `xgb_model` / `xgb_cutoff` — documented in
[Consensus model](#consensus-model).

**Per-caller tuning** blocks exist for CoNIFER (`conifer_*`), XHMM (`xhmm_*`),
jabCoNtool (`jabCoNtool_*`) and Control-FREEC (`freec_*`).

## What the pipeline does

1. **Preprocessing**
   - `BED_PREPARATION` — binned genome + GC profile from the FASTA/dict
   - `BED_GTF_ANNOTATION` — annotate the panel BED with gene/biotype from the GTF
   - `PER_REGION_COVERAGE_CALC` — per-region coverage per sample (used in the final tables)
2. **Callers** (each gated by its toggle, each writing one per-sample calls file)

   | Caller | Subworkflow | Calls file consumed downstream |
   |---|---|---|
   | CNVkit | `cnvkit` | `<sample>_cnvkit.tsv` |
   | GATK gCNV | `gatk` | `<sample>_gatk.vcf.gz` (segments) |
   | ExomeDepth | `ExomeDepth` | `<sample>_ExomeDepth.tsv` |
   | panelcn.MOPS | `panelcnMOPS` | `<sample>_panelcnMOPS.tsv` |
   | cn.MOPS | `cnMOPS` | `<sample>_cnMOPS.tsv` |
   | XHMM | `XHMM` | `<sample>_xhmm.tsv` |
   | CoNIFER | `conifer` | `<sample>_conifer.tsv` |
   | Control-FREEC | `controlFREEC` | `<sample>_freec.tsv` |
   | ECOLE (deep learning) | `ECOLE` | per-sample text output, published only (not combined, see [Pitfalls](#pitfalls-and-caveats)) |
   | jabCoNtool | `jabcontool` | normalised per-sample varcalls (added last) |

   Calls are accumulated per sample with `Utils.mixAndCollectVarcalls()` into
   `ch_all_varcalls` → `ch_all_varcalls_for_merging` = `[meta, [all caller files]]`.
3. **Combination** — either `MERGE_VARIANT_CALLS` (R merge & smooth) or
   `XGBOOST_CONSENSUS_SCORE` (score with the trained model, filter by `xgb_cutoff`, refine
   the survivors to supported panel exons).
4. **Annotation** — `CLASSIFY_AND_ANNOTATE` (ClassifyCNV + `cnvAnnotateCNVkit.R`).
5. **Reporting** — `FINAL_FORMATTING_TABLES` produces the cohort tables and HTML report.
6. **Optional** — `PSEUDOGENE_ANALYSIS`, gated by `run_pseudogene` and driven by
   `gene_pseudogene_bed`.

## Output layout

Relative to the launch directory:

```
structural_varcalls/
  <sample>/                       per-caller raw calls
    cnvkit/ exomeDepth/ panelcnMOPS/ gatk/ XHMM/ conifer/ freec/ ...
  all_samples/                    cohort-level caller outputs (ECOLE, jabCoNtool)
  bed_annotation/                 annotated panel BED
  cohort_results/
    all_samples_merged.tsv
    all_samples_annotated.tsv
    all_samples_combined.tsv
    all_samples_smoothed.tsv
    final_report.html
    coverage/*.tsv
coverage_tracks/<sample>/         per-region coverage
cohort_data/                      reusable cohort references (cn.MOPS etc.)
pseudogene/final_classification/  optional pseudogene analysis
```

Per-sample combined calls (`merged_variants/*`) and annotated calls
(`classified_and_annotated_CNVs/*`) are not published directly; they flow into
`cohort_results/`.

---

# Consensus model

`modules/xgboost_consensus_model/` replaces the rule-based merge with a learned model that
assigns a **probability that a candidate region is a real CNV**, given which callers detected
it, what type it is, and how wide it is.

Everything belonging to the model — the deployed scorer *and* the offline research/training
code that produced it — lives under that one directory:

```
modules/xgboost_consensus_model/
  cnv_consensus_model.py                  shared library: per-caller parsers, overlap logic,
                                          candidate merging, ground-truth loading
  score_samples_xgb.py                    scoring entry point, one sample per invocation
  refine_regions.py                       breaks scored candidates down to supported panel exons
  xgboost_tuned_no_freec_cnmops.json      the model shipped with the repo (XGBoost native format)
  main.nf                                 XGBOOST_CONSENSUS_SCORE
  env.yaml                                python 3.12 + numpy/pandas/xgboost
  modelling/                              offline training + evaluation (not run by the pipeline)
    merge_from_bronco.py                  build the training table from the BRONCO cohort
    src/                                  load → features → baseline/eda → tune → train
    models/                               trained model artifacts
    reports/metrics/                      metrics for every training phase
    xgboost_pipeline/                     standalone pipeline: score pre-computed caller output
```

**Paths in the rest of this section are relative to `modules/xgboost_consensus_model/`.**

The model is **trained offline**, not by the pipeline — see [Training the model](#training-the-model).

## Candidate regions

Calls from all callers are clustered by chromosome and variant type (DEL/DUP) into candidate
regions. A call counts as overlapping when the overlap covers at least `MIN_OVERLAP_FRAC = 0.5`
of the *smaller* interval, so a gene-level segment containing a single exon still matches it.

**The model never chooses coordinates.** `_merge_candidates()` does, and it takes the **union**
of the cluster: `start = min(start)`, `end = max(end)` over its members. A region is therefore
as wide as the widest caller that supports it. These union coordinates are what the model
scores; they are only narrowed afterwards, by `refine_regions.py`, for the candidates that pass
the cutoff (see [Scoring](#scoring-use_consensus--true), step 5).

To stop one caller inflating everything, calls larger than `MAX_CLUSTER_SIZE_RATIO = 50` × the
median call size in a cluster are treated as **coordinate outliers**: they still count toward
`n_callers` and `callers_supporting`, but their boundaries are excluded. That way a caller
reporting a gene-level deletion can vote for a candidate without stretching it to megabase scale.

Clustering is greedy and seeded in sorted-by-start order, and each absorbed call widens the
interval that subsequent calls are tested against — so clusters chain: if A–B overlap and B–C
overlap, C joins even when A and C do not.

## Features and model

Ten columns, in this exact order (`FEATURE_COLUMNS` in `score_samples_xgb.py`, which must stay
in sync with `modelling/src/train_no_freec_cnmops.py`):

| Feature | Meaning |
|---|---|
| `cnvkit`, `exomeDepth`, `panelcnMOPS`, `XHMM`, `conifer`, `gatk` | binary — did that caller support the candidate |
| `agreement_count` | number of supporting callers, counted over **all 8** of `ALL_CALLERS` |
| `type_DUP` | `1` for DUP, `0` for DEL |
| `log_span_width` | `log1p(end - start)` |
| `type_conflict` | do DEL and DUP both exist at the same `(chr, start, end)` |

`freec` and `cnMOPS` are dropped as model inputs but **still counted in `agreement_count`** —
that is how the training data was built, and "simplifying" it to 6 would silently mismatch the
model. Per-caller confidence scores are not used; only the detection pattern, the variant type,
and the span.

**Model.** `xgboost.XGBClassifier`, 500 trees, `max_depth=4`, `learning_rate=0.05`,
`scale_pos_weight≈7.56`, loaded from XGBoost native JSON. Out-of-fold row precision `0.884` /
recall `0.928` / F2 `0.919`, versus `0.784` / `0.922` for the "gatk fired" baseline — see
[`modelling/reports/metrics/no_freec_cnmops.md`](modules/xgboost_consensus_model/modelling/reports/metrics/no_freec_cnmops.md).
Its own ≥0.90-precision operating threshold is **`0.7499`**.

**Callers.** The parsers understand
`cnvkit, freec, exomeDepth, panelcnMOPS, cnMOPS, XHMM, conifer, gatk` (`ALL_CALLERS`).
Files are located by pattern:

| Caller | Pattern |
|---|---|
| cnvkit | `*_cnvkit.tsv` |
| freec | `*_freec.tsv` |
| exomeDepth | `*_ExomeDepth.tsv` |
| panelcnMOPS | `*_panelcnMOPS.tsv` |
| cnMOPS | `*_cnMOPS.tsv` |
| XHMM | `*_xhmm.tsv` |
| conifer | `*_conifer.tsv` |
| gatk | `*_gatk.vcf.gz` (falls back to `*_gatk.vcf`) |

Both directory layouts are recognised — `<sample>/<caller>/<file>` and the flat
`<sample>/<sample>_<caller>.tsv` that this pipeline produces.

## Scoring (`use_consensus = true`)

```groovy
use_consensus = true
xgb_cutoff    = 0.7499   // or null to keep every scored candidate
xgb_model     = null     // or '/path/to/your_model.json'
```

| Param | Default | Meaning |
|---|---|---|
| `use_consensus` | `true` | `true` = score with the model instead of `MERGE_VARIANT_CALLS` |
| `xgb_model` | `null` | pre-trained XGBoost `.json`; `null` = the one shipped in the module |
| `xgb_cutoff` | `0.7499` | score threshold; the model's ≥0.90-precision operating point (precision `0.9002` / recall `0.9191`). `null` keeps every scored candidate |
| `smooth_gap_bp` | `10000` | gap for re-merging refined exons (step 5); also used by `FINAL_FORMATTING_TABLES` |

`min_number_callers` and `size_threshold` belong to the rule-based path and are **ignored**
when `use_consensus = true` — `xgb_cutoff` is the only filter.

### What happens, per sample

`XGBOOST_CONSENSUS_SCORE` runs once per sample, in place of `MERGE_VARIANT_CALLS`. Everything
upstream (the callers) and downstream (`CLASSIFY_AND_ANNOTATE`, `FINAL_FORMATTING_TABLES`) is
unchanged.

1. **Stage the calls.** The sample's caller files are linked into `sv_dir/<sample>/`, and the
   calls of all 8 callers in `ALL_CALLERS` are parsed.
2. **Build candidates.** `_merge_candidates()` clusters the calls into union-coordinate
   candidate regions, as described in [Candidate regions](#candidate-regions).
3. **Score.** `score_samples_xgb.py` builds the ten [features](#features-and-model) for each
   candidate and writes the model's probability as `consensus_score`.
4. **Filter.** Candidates with `consensus_score >= xgb_cutoff` continue; the rest are dropped
   from the annotated and final tables. The complete scored set is still kept in
   `_all_scored_candidates.tsv`.
5. **Refine to supported exons.** `refine_regions.py` intersects each surviving candidate with
   the panel BED and keeps only the exons that a member call actually overlaps. A member call
   doesn't count as exon-level support when it is a size outlier (`MAX_CLUSTER_SIZE_RATIO`) or
   spans exons of more than `MAX_SUPPORT_GENES = 2` panel genes. Without this step, a single
   multi-megabase segment would make every exon it crosses an affected exon. Exons kept for
   the same sample, chromosome and type are re-merged across gaps of up to `smooth_gap_bp`.
   A candidate with no panel overlap, or no exon-level support, is passed through **unchanged**
   rather than dropped. Losing a called region silently would be worse than reporting a coarse
   one.
6. **Reshape.** The refined regions are rewritten to the exact 13-column `MERGE_VARIANT_CALLS`
   schema as `_smoothed_variants.tsv/.bed`. This is what `CLASSIFY_AND_ANNOTATE` receives, so
   genes, classification and score are computed on the refined coordinates, not the union
   span.

### Outputs

Written into `merged_variants/` in the task's work directory (not published; they flow into
`cohort_results/`):

| File | Content |
|---|---|
| `<sample>_all_scored_candidates.tsv` | **every** scored candidate, unfiltered and unrefined. Use this to re-threshold or audit |
| `<sample>_refined_regions.tsv` | the final region set after steps 4–5 |
| `<sample>_smoothed_variants.tsv` | the refined regions in the 13-column `MERGE_VARIANT_CALLS` schema, passed to annotation |
| `<sample>_smoothed_variants.bed` | 1-based START, so `VariantID`s match ClassifyCNV's |
| `<sample>_merged_target_consensus.tsv/.bed` | the cutoff-filtered candidates at union coordinates (before refinement) |
| `<sample>_per_exon_matrix.tsv` | copy of `_all_scored_candidates.tsv`, for downstream compatibility |

Raw model output columns are
`sample | chr | start | end | type | n_callers | callers_supporting | consensus_score | cn_label`.
With `xgb_cutoff` set, `score_samples_xgb.py` also writes
`*_smoothed_variants_filtered_<cutoff>.tsv`, and steps 5–6 read from that file. When the
cutoff is set and candidates exist but the file is missing, the task fails rather than passing
on unfiltered calls labelled as filtered.

The same module backs the standalone
[`modelling/xgboost_pipeline/`](modules/xgboost_consensus_model/modelling/xgboost_pipeline/README.md), which scores
pre-computed caller output without re-running the callers.

## Training the model

Training is **offline only** — there is no in-pipeline training process. The pipeline scores
with a pre-trained model and nothing else. The full path, all under `modelling/`:

| Step | Script | Produces |
|---|---|---|
| 1. Merge caller calls across the training cohort | `merge_from_bronco.py` | `all_samples_merged.tsv` |
| 2. Build the feature matrix | `src/features.py` | feature table + grouped CV folds |
| 3. Baseline / EDA (optional) | `src/baseline.py`, `src/eda.py` | `reports/metrics/baseline.*` |
| 4. Tune hyperparameters | `src/tune.py` | `reports/metrics/tuning.*` |
| 5. Fit the shipped model | `src/train_no_freec_cnmops.py` | `models/xgboost_tuned_no_freec_cnmops.json` + `reports/metrics/no_freec_cnmops.*` |

Splits are **grouped by sample** (`StratifiedGroupKFold`), so no sample appears in both train
and test. Metrics are reported out-of-fold, and `F2` is the selection criterion throughout —
recall is weighted over precision, on the view that a missed CNV costs more than one to review.

To deploy a newly trained model, copy it from `modelling/models/` up to the module root as
`xgboost_tuned_no_freec_cnmops.json` (what `xgb_model = null` uses), or point `xgb_model` at
it. **If you change `FEATURE_COLUMNS`, change it in both
`train_no_freec_cnmops.py` and `score_samples_xgb.py`** — the column order must match the
model's stored `feature_names` or scoring silently mispredicts.

## Pitfalls and caveats

**Ground-truth sample names must match exactly.** Ground truth is keyed by the `sample` *column
inside* the TSV, not by the filename, and is intersected with sample directory names — which in
this pipeline are `meta.sample_name`, i.e. the full BAM basename
(`007_CZE3xBRCA1328_run137.picard.sorted.RG.rmdup`), not a short label. A mismatch yields
0 usable samples.

**Negatives assume complete ground truth.** Every caller call that does not overlap a
ground-truth variant becomes a negative training example. That is correct for a simulated cohort
where the injected variants are the whole truth, but if the ground truth lists only a handful of
known variants per sample, every other genuine CNV is labelled a false positive and the model
learns to distrust the callers.

**Coordinates are only as good as the callers.** The model scores union coordinates, so a single
over-called segment widens the region it sees (`log_span_width`). `agreement_count` also counts
a caller that agrees on the *event* even when it disagrees badly on the *boundaries*.
`refine_regions.py` narrows the reported regions to supported exons. It does not change what
was scored, and candidates outside the panel BED keep their union span. If breakpoint precision
matters, `MIN_OVERLAP_FRAC`, `MAX_CLUSTER_SIZE_RATIO` and `MAX_SUPPORT_GENES` are the knobs.
Retraining will not help, since the model never sees a coordinate choice to learn from.

**Keep all 8 callers enabled.** The model's direct features are six callers, but `freec` and
`cnMOPS` still feed it. Both count in `agreement_count`, and both take part in candidate
merging, which decides which candidates exist, their boundaries (`log_span_width`) and
`type_conflict`. Turning any of the eight off shifts every candidate's features away from what
the model saw in training. The `0.7499` cutoff then no longer gives the ~0.90 precision the
report states.

**A failed scoring step fails the task.** A sample with no caller calls gets a header-only
table from `score_samples_xgb.py` and carries on with zero CNVs. Any non-zero exit is therefore
a real error, and it stops the task rather than being reported as a sample with no CNVs:
a conda env without `xgboost` ≥ 3.4, an unreadable model, a parser error. The task's
`.command.err` has the traceback (`nextflow log <run_name> -f name,workdir`).

**ECOLE is not combined.** ECOLE calls the whole cohort in one task and its output carries no
per-sample label, so it is kept out of `ch_all_varcalls`. Neither combiner reads it. With
`use_ecole = true` its calls are only published to `structural_varcalls/all_samples/`.

**`XGBOOST_CONSENSUS_SCORE` is `cache false`.** Scoring re-runs on every `-resume`.

---

## Cohort data reuse

cn.MOPS, panelcn.MOPS and ExomeDepth build a cohort reference from all samples in the run. To
reuse a pre-built one instead:

```groovy
use_panelcnmops_cohortdata = true
cohort_panelcnmops_ref     = '/path/HyperExomeV2_GRCh38_panelcnmops.RData'
use_cnmops_cohortdata      = true
cohort_cnmops_ref          = '/path/HyperExomeV2_GRCh38_cnmops.RData'
use_exomedepth_cohortdata  = true
cohort_exomedepth_ref      = '/path/HyperExomeV2_GRCh38_exomedepth.RData'
```

With these `false`, the run's own BAMs are the cohort — which needs enough samples for the
statistics to be meaningful.

## Resources and executor tuning

`nextflow.config` sets a local executor (`cpus = 25`, `memory = '90 GB'`, `queueSize = 2`) and
per-process `maxForks`/`memory` in the `process { withName: ... }` block. The cohort-preparation
and coverage steps are the memory-hungry ones (`90 GB`, `maxForks = 1`); the per-sample callers
run with `maxForks = 4`. Adjust both the `executor` block and the per-process overrides together
— raising `queueSize` without lowering the per-process memory will over-commit the machine.

`local.config` holds site-specific overrides (tmp directories, `conda.cacheDir`, `scratch`) and
is applied with `-c local.config`.

