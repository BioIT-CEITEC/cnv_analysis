# cnv_analysis — germline CNV/SV calling pipeline

A Nextflow (DSL2) pipeline that runs a panel of independent CNV callers over targeted
(panel/WES) BAM files, combines their calls into a per-sample consensus, annotates and
classifies the result, and emits cohort-level tables plus an HTML report.

Two combination strategies are available:

| Strategy | Enabled by | Implementation |
|---|---|---|
| Rule-based merge & smooth (default) | `use_consensus = false` | `bin/merging_and_smoothing.R` via `MERGE_VARIANT_CALLS` |
| **Trained consensus model** | `use_consensus = true` | `modules/consensus_model/` — logistic model with per-caller weights, **trainable from ground-truth data** (see [Consensus model](#consensus-model)) |

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
  - [How the model works](#how-the-model-works)
  - [Scoring mode](#scoring-mode)
  - [Training mode](#training-mode)
  - [Ground-truth format](#ground-truth-format)
  - [Training parameters](#training-parameters)
  - [Training outputs](#training-outputs)
  - [Reusing a trained model](#reusing-a-trained-model)
  - [Pitfalls and caveats](#pitfalls-and-caveats)
- [Cohort data reuse](#cohort-data-reuse)
- [Resources and executor tuning](#resources-and-executor-tuning)
- [Known limitations](#known-limitations)

---

## Requirements

- **Nextflow** ≥ 22.10 — DSL2, and the consensus-model wiring uses `error` and directory
  `stageAs` patterns
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
  WGS.nf                    WGS variant (currently not included from main.nf)
subworkflows/<tool>/main.nf per-caller orchestration
modules/<tool>/main.nf      individual processes (+ env.yaml, helper scripts)
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

**Combination**: `use_consensus`, plus the whole `consensus_*` block — documented in
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
   | ECOLE (deep learning) | `ECOLE` | per-sample text output |
   | jabCoNtool | `jabcontool` | normalised per-sample varcalls (added last) |

   Calls are accumulated per sample with `Utils.mixAndCollectVarcalls()` into
   `ch_all_varcalls` → `ch_all_varcalls_for_merging` = `[meta, [all caller files]]`.
3. **Combination** — either `MERGE_VARIANT_CALLS` (R merge & smooth) or
   `CONSENSUS_MODEL_SCORE` (trained model).
4. **Annotation** — `CLASSIFY_AND_ANNOTATE` (ClassifyCNV + `cnvAnnotateCNVkit.R`).
5. **Reporting** — `FINAL_FORMATTING_TABLES` produces the cohort tables and HTML report.
6. **Optional** — `PSEUDOGENE_ANALYSIS` (`run_pseudogene`, see
   [Known limitations](#known-limitations)).

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
consensus_model/                  trained consensus model + training statistics
pseudogene/final_classification/  optional pseudogene analysis
```

Per-sample combined calls (`merged_variants/*`) and annotated calls
(`classified_and_annotated_CNVs/*`) are not published directly; they flow into
`cohort_results/`.

---

# Consensus model

`modules/consensus_model/` replaces the rule-based merge with a learned model that assigns a
**probability that a candidate region is a real CNV**, given which callers detected it and how
confident each one was.

```
modules/consensus_model/
  cnv_consensus_model.py      shared library: parsers, overlap logic, candidate merging,
                              BoundedConsensusModel, feature matrix
  train_consensus_model.py    training entry point (ground truth → model)
  score_samples.py            scoring entry point (model → per-sample scored calls)
  cnv_consensus_model.pkl     the model shipped with the repo (used when not training)
  cnv_consensus_model.json    its weights/metrics in readable form
  cnv_train_*.tsv             training statistics from the shipped model
  main.nf                     CONSENSUS_PREPARE_TRAINING_DIR, CONSENSUS_MODEL_TRAIN,
                              CONSENSUS_MODEL_SCORE
  env.yaml                    python 3.12 + numpy/pandas/scikit-learn/scipy
```

## How the model works

**Candidate regions.** Calls from all callers are clustered by chromosome and variant type
(DEL/DUP) into candidate regions. A call counts as overlapping when the overlap covers at least
`MIN_OVERLAP_FRAC = 0.5` of the *smaller* interval, so a gene-level segment containing a single
exon still matches it. Calls larger than `MAX_CLUSTER_SIZE_RATIO = 50` × the median call size in
a cluster still count toward `n_callers` but do not extend the candidate's boundaries — this
stops one caller's multi-megabase segment from inflating an exon-level candidate.

**Features.** Two per caller: `[hit, score]`, where `hit ∈ {0,1}` is detection of the candidate
and `score` is that caller's own confidence for its best matching call (0 when there is no hit).

**Model.** `BoundedConsensusModel` = `StandardScaler` + logistic regression, with the constraint
that **every coefficient is non-negative and capped**. A caller can therefore only ever add
evidence, never subtract it, and no single caller can dominate:

```
consensus_score(region) = sigmoid( intercept + Σ_c  w_c · features_c(region) )
```

Two capping modes:

- `per_caller` (default) — `hit` + `score` coefficients of one caller sum to ≤ `weight_cap`
- `per_coef` — each coefficient individually ≤ `weight_cap` (a caller can reach 2 × cap)

**Callers.** The parsers understand
`cnvkit, freec, exomeDepth, panelcnMOPS, cnMOPS, XHMM, conifer, gatk`
(`ALL_CALLERS`). The default trained subset is
`exomeDepth, panelcnMOPS, XHMM, conifer, gatk` (`DEFAULT_CALLERS`). Files are located by pattern:

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

## Scoring mode

```groovy
use_consensus    = true
consensus_cutoff = null     // or e.g. 0.5
consensus_model  = null     // or '/path/to/your_model.pkl'
```

`CONSENSUS_MODEL_SCORE` runs once per sample. It stages that sample's caller files into
`sv_dir/<sample>/`, runs `score_samples.py`, and writes into `merged_variants/`:

| File | Content |
|---|---|
| `<sample>_smoothed_variants.tsv` | all scored candidates, reshaped to the exact `MERGE_VARIANT_CALLS` schema (13 columns) so the annotation step is agnostic to which combiner ran |
| `<sample>_smoothed_variants.bed` | 1-based START, so `VariantID`s match ClassifyCNV's |
| `<sample>_merged_target_consensus.tsv/.bed` | the cutoff-filtered set (or all calls when no cutoff) |
| `<sample>_per_exon_matrix.tsv` | copy of the scored file, for downstream compatibility |

Raw model output columns are
`sample | chr | start | end | type | n_callers | callers_supporting | consensus_score | cn_label`.
`consensus_cutoff` additionally makes `score_samples.py` write
`*_smoothed_variants_filtered_<cutoff>.tsv`, which becomes the source of
`_merged_target_consensus.*`.

The process also restores the callers' original chromosome naming: `norm_chr()` always
normalises internally to `chrN`, so the module detects whether the input files used a `chr`
prefix and strips it again if they did.

## Training mode

Set `consensus_train = true` to fit a new model from ground-truth data. Training runs **before**
scoring, and when `use_consensus = true` the freshly trained model is the one used to score the
run. Training also works with `use_consensus = false` — the model is trained and published, and
the run still combines calls with the rule-based merge.

Two processes are involved:

- **`CONSENSUS_PREPARE_TRAINING_DIR`** — only when training on this run's own calls. Copies each
  sample's caller files into a directory named after the sample, which is the layout the trainer
  expects.
- **`CONSENSUS_MODEL_TRAIN`** — stages all sample directories under `sv_dir/`, runs
  `train_consensus_model.py`, and publishes the model to `consensus_model/`.

### Two training cohorts

**A. External cohort (recommended).** Point the pipeline at a directory of previously called
samples — typically a simulated cohort where the injected variants *are* the complete truth:

```groovy
use_consensus          = true
consensus_train        = true
consensus_train_sv_dir = '/data/model_training/structural_varcalls_simulated_BR'
consensus_train_gt_dir = '/data/model_training/selected_regions_simulated'
```

Each immediate sub-directory of `consensus_train_sv_dir` is one training sample; the names
`all_samples`, `selected_regions`, `merged_variants` and `selected_regions_simulated` are
skipped automatically, so pointing at an existing `structural_varcalls/` directory works.

**B. This run's own calls.** Leave `consensus_train_sv_dir` unset and the samples being analysed
become the training cohort:

```groovy
use_consensus          = true
consensus_train        = true
consensus_train_gt_dir = '/data/truth/my_cohort'
```

Only sensible when the ground truth for those samples is (near-)complete — see
[Pitfalls and caveats](#pitfalls-and-caveats).

### What training does

1. Load ground truth; merge exon-level entries into variant-level events (same chromosome and
   type, gap ≤ 10 kb).
2. Intersect labelled samples with the sample directories actually present.
3. Detect which callers are present; restrict to `consensus_train_callers` or to the
   `DEFAULT_CALLERS` that exist in the data.
4. Split at **sample** level, stratified by "has any ground-truth variant"
   (`consensus_train_test_size`, `consensus_train_seed`). The test split is recorded but not
   consumed by training — it exists so held-out samples can be evaluated separately.
5. Per-caller TP/FP/FN/precision/recall/F1 on the training samples, overall and split by DEL/DUP.
6. Build the feature matrix — **positives** = one row per ground-truth variant; **negatives** =
   each caller call that overlaps no ground-truth variant. Training aborts below 5 positives.
7. Fit with 5-fold stratified CV (ROC-AUC), then refit on all training rows.
8. Derive the **max-F1 threshold** on the training data as `suggested_cutoff`.

## Ground-truth format

One TSV per sample in `consensus_train_gt_dir`, named `<SAMPLE>_selected_regions.tsv`, with a
header row:

| Column | Meaning |
|---|---|
| `sample` | sample identifier — **must equal the sample directory name** |
| `chrom` | chromosome (with or without `chr`; normalised internally) |
| `start` | 1-based start |
| `end` | end |
| `name` | free text (gene/exon label), not used for matching |
| `variant_type` | free text |
| `type` | **`DEL` or `DUP`** — anything else is ignored |

Rows whose `type` is not DEL/DUP are dropped silently. Empty files and dot-files are skipped.
Exon-level rows belonging to one event should simply be listed individually; the 10 kb merge
step reassembles them.

Example:

```
sample	chrom	start	end	name	variant_type	type
BR-2604	chr17	41243452	41246877	BRCA1_ex11	simulated	DEL
BR-2604	chr17	41247862	41249306	BRCA1_ex10	simulated	DEL
BR-2604	chr13	32900238	32900420	BRCA2_ex3	simulated	DUP
```

## Training parameters

| Param | Default | Meaning |
|---|---|---|
| `consensus_train` | `false` | master switch for training |
| `consensus_train_gt_dir` | `null` | **required** when training — ground-truth directory |
| `consensus_train_sv_dir` | `null` | external training cohort; `null` = use this run's calls |
| `consensus_train_callers` | `null` | space-separated subset, e.g. `'exomeDepth panelcnMOPS XHMM conifer gatk'` |
| `consensus_train_weight_cap` | `1.0` | upper bound of each non-negative caller weight |
| `consensus_train_cap_mode` | `'per_caller'` | `per_caller` or `per_coef` (see above) |
| `consensus_train_l2` | `0.01` | L2 penalty on the coefficients |
| `consensus_train_test_size` | `0.2` | fraction of samples held out (`0.0` = train on everything) |
| `consensus_train_seed` | `42` | split seed |
| `consensus_cutoff` | `null` | scoring threshold; when set it is also recorded in the trained model in place of the max-F1 cutoff |
| `consensus_model` | `null` | pre-trained `.pkl` to score with; ignored while `consensus_train = true` |

## Training outputs

Published to `consensus_model/`:

| File | Content |
|---|---|
| `cnv_consensus_model.pkl` | the pickled model — what scoring loads |
| `cnv_consensus_model.json` | caller weights, per-feature weights, intercept, cutoff, callers, CV/train metrics, and the train/test sample lists |
| `cnv_train_caller_stats.tsv` | per-caller TP/FP/FN/precision/recall/F1 (training samples) |
| `cnv_train_caller_stats_by_type.tsv` | the same, split into DEL and DUP |
| `cnv_train_pr_curve.tsv` | consensus precision/recall/F1 across thresholds |
| `cnv_train_train_test_split.tsv` | which sample went to train vs test |

The training log prints the usable-sample count, the callers actually used, the per-caller
table, the metrics and the final weights — **read it before trusting the model**. Key lines:

```
  N samples usable (present in GT + simulated)
  Using callers: exomeDepth, panelcnMOPS, XHMM, conifer, gatk
  Shape: (rows, features)  |  positives: P  negatives: N
  Suggested cutoff (max-F1 on train): 0.XXXX
```

Judge the result from `cv_roc_auc_mean` (generalisation) versus `train_roc_auc` (fit), and from
`cnv_train_pr_curve.tsv` when choosing a `consensus_cutoff`.

## Reusing a trained model

Training every run is wasteful. Once you are happy with a model, either point the pipeline at
the published copy:

```groovy
consensus_train = false
use_consensus   = true
consensus_model = '/path/to/consensus_model/cnv_consensus_model.pkl'
```

or promote it to the repo default by copying it (and its `.json`) over
`modules/consensus_model/cnv_consensus_model.pkl`, which is what `consensus_model = null` uses.

The `.pkl` is a Python pickle of `BoundedConsensusModel`, so it must be loaded with a compatible
`cnv_consensus_model.py` and a scikit-learn close to the one that produced it. Keep the `.json`
alongside it — it is the human-readable record of how that model was trained.

## Pitfalls and caveats

**Sample names must match exactly.** Ground truth is keyed by the `sample` *column inside* the
TSV, not by the filename, and is intersected with the sample directory names — which in this
pipeline are `meta.sample_name`, i.e. the full BAM basename
(`007_CZE3xBRCA1328_run137.picard.sorted.RG.rmdup`), not a short label. A mismatch yields
0 usable samples and a confusing failure inside the train/test split rather than a clear error.
Check the `N samples usable` line first whenever training fails.

**Negatives assume complete ground truth.** Every caller call that does not overlap a
ground-truth variant becomes a negative training example. That is correct for a simulated cohort
where the injected variants are the whole truth, but if the ground truth lists only a handful of
known variants per sample, every other genuine CNV is labelled a false positive and the model
learns to distrust the callers. For real cohorts, prefer an external simulated/curated training
set (`consensus_train_sv_dir`).

**Training on the run you score is optimistic.** Cohort B trains on the same calls it then
scores; the resulting scores are not an unbiased estimate of performance. Use the held-out split
or a separate cohort for any performance claim.

**Callers present ≠ callers expected.** The model is only trained on callers whose files exist
in the training cohort, and scoring uses the caller list stored in the model. If you train with
a cohort produced by a different toggle set than the run you later score, the weights will not
cover the callers you actually have. Compare `Using callers:` in the training log with your
`use_*` toggles.

**`CONSENSUS_MODEL_SCORE` is `cache false`.** Scoring re-runs on every `-resume`; training does
not.

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

## Known limitations

- **Pseudogene branch is not wired.** `params.run_pseudogene` and `params.organism_gene_bed` are
  not declared, and `PSEUDOGENE_ANALYSIS` is called with an undefined `ch_organism_gene_bed`.
  Enabling it as-is will fail.
- **WGS workflow is disabled.** `workflows/WGS.nf` exists but its include is commented out in
  `main.nf`; only `PANEL_WES` runs.
- **ECOLE fine-tuning is not exposed.** `modules/ECOLE/ECOLE-0.2/` ships
  `ECOLE_finetune.py`, `finetune_preprocess_sample.py` and `create_dataset.py`, but only
  read-depth preprocessing and calling are wired into the workflow; ECOLE always uses its
  pre-trained weights. The scripts also hardcode paths
  (`./finetune_example_data/ground_truth_labels/`, `./processed_finetuning_dataset`) that would
  need parameterising first.
- **Unused modules.** `delly`, `purple`, `vardict`, `sage`, `svprep`, `cobalt`, `amber` and
  `modules/backup/*` are present but not part of `PANEL_WES`.
- **Heavy ECOLE environment.** Pinned to Python 3.8 / CUDA 10.1 via a full Anaconda export;
  expect a long first-run environment build. `use_ecole = false` skips it.
