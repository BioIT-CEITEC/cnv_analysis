# XGBoost Consensus Scoring Pipeline

Standalone Nextflow pipeline that scores pre-computed CNV caller output with
the tuned XGBoost consensus model (no `freec`/`cnMOPS`), then annotates and
produces cohort-level tables + an HTML report.

**Does not run the callers.** Point it at caller output that already exists.

## Input layout

One directory per sample, holding that sample's caller output files:

```
<sample_dir>/
  cnvkit/...        (*_cnvkit.tsv)
  freec/...         (*_freec.tsv)
  exomeDepth/...    (*_ExomeDepth.tsv)
  panelcnMOPS/...   (*_panelcnMOPS.tsv)
  cnMOPS/...        (*_cnMOPS.tsv)
  XHMM/...          (*_xhmm.tsv)
  conifer/...       (*_conifer.tsv)
  gatk/...          (*_gatk.vcf.gz)
```

Caller-named subdirectories (above) or flat `<sample>_<caller>.*` files both
work. Not all 8 callers are required — whichever are present get used.

## Point it at your data

Add each sample to `params.samples` — either directly in
`nextflow.config`, or (**required for real/private sample data — see
below**) in a separate `-c` override file:

```groovy
params {
    samples = [
        'my_sample': [sample_name: 'my_sample', sv_dir: '/path/to/my_sample_caller_dir'],
    ]
}
```

`sv_dir` defaults to `sv_dir/<sample_name>/` next to this pipeline if
omitted.

**Real/private sample data:** never put real sample names or paths in
`nextflow.config` (it's tracked in git). Put them in a separate config file
instead and gitignore it, e.g. `test_prague.config`:

```groovy
params {
    samples = [
        'sample1': [sample_name: 'sample1', sv_dir: '/real/path/to/sample1'],
    ]
}
```

Then run with an extra `-c` (see [Run it](#run-it) for `$PIPE`):
```bash
nextflow run $PIPE/main.nf -c local.config -c $PIPE/test_prague.config -resume
```

## Run it

```bash
cd <repo root>
PIPE=modules/xgboost_consensus_model/modelling/xgboost_pipeline

# 1. Dry run first -- checks wiring, no real work
nextflow run $PIPE/main.nf -c local.config -stub -resume

# 2. Real run
nextflow run $PIPE/main.nf -c local.config -resume
```

## Config

| Param | Default | Meaning |
|---|---|---|
| `samples` | `[:]` | See above (defined here) |
| `xgb_model` | `null` | Override the model path; default is `<repo root>/modules/xgboost_consensus_model/xgboost_tuned_no_freec_cnmops.json` |
| `xgb_cutoff` | `0.7499` | Score threshold to filter on: the model's own ≥0.90-precision threshold, inherited from the repo-root config. To keep every scored candidate, set `params.xgb_cutoff = null` in a `-c` override file. Don't pass `--xgb_cutoff null` on the command line, which may arrive as the string `"null"`. The unfiltered set is always in `*_all_scored_candidates.tsv` |

`xgb_model` / `xgb_cutoff` are declared in the repo-root `nextflow.config`.
This pipeline lives inside `modules/xgboost_consensus_model/`, the same module
the main `CNV_ANALYSIS` pipeline scores with (`use_consensus = true`), so both
entry points run the identical process and read the same two params.

All `organism_*` reference params (genome, GTF, DNA panel, ...) are
inherited from the repo-root `nextflow.config` — nothing to set here unless
overriding. `organism_dna_panel` is **required**: scoring breaks each
candidate down to the supported exons of that panel BED (`refine_regions.py`)
before annotation. If it is unset, the panel channel is empty and no sample
is ever scored. The panel must be the one the caller output was produced
against.

## Output

| What | Where |
|---|---|
| Cohort tables + HTML report | `structural_varcalls/cohort_results/` (relative to where you ran `nextflow run` from) |
| Per-sample scored candidates / annotated calls | Nextflow `work/` only — no `publishDir` on those steps (same as the main pipeline's `XGBOOST_CONSENSUS_SCORE`/`CLASSIFY_AND_ANNOTATE`). Find a specific task with `nextflow log <run_name> -f hash,name,workdir`. |
