# AGENT.md — working context for AI agents

Read this before changing anything. [README.md](README.md) is the user-facing reference;
this file covers what an agent needs to edit the repo safely. Where they overlap, the README has
the detail.

## What this repo is

A Nextflow DSL2 pipeline for **germline CNV calling on targeted (panel/WES) BAMs**. It does four
things in order:

1. Runs up to 10 independent CNV callers.
2. Combines their calls per sample, either with the rule-based `MERGE_VARIANT_CALLS` or with
   the XGBoost consensus model `XGBOOST_CONSENSUS_SCORE` (`use_consensus = true`, the default).
3. Annotates the combined calls with ClassifyCNV.
4. Writes cohort tables and an HTML report.

Call chain: `main.nf` → `CNV_ANALYSIS` → `PANEL_WES` in `workflows/PANEL.nf`. That runs
`subworkflows/<caller>/` → `modules/<tool>/` → helper scripts in `bin/`.

The XGBoost model is **trained offline**. The pipeline only scores. Everything model-related
lives in `modules/xgboost_consensus_model/`:

- the deployed scorer: `score_samples_xgb.py`, `refine_regions.py`, `cnv_consensus_model.py`,
  `main.nf`
- the shipped model: `xgboost_tuned_no_freec_cnmops.json`
- the research and training code under `modelling/`: `src/` and `merge_from_bronco.py`
- a standalone scoring pipeline: `modelling/xgboost_pipeline/`

## Environment facts

- **The pipeline runs on the CEITEC cluster, not on the developer's Mac.** Reference,
  sample and tmp paths in `nextflow.config` are `/mnt/data/ceitec_cfg2/...` cluster paths, and
  Nextflow is not installed locally. You cannot run the pipeline here. Say so instead of
  claiming a change works.
- **Local Python:** `~/miniforge3/envs/cnv_consensus_ml/bin/python` has pandas but **not
  xgboost**. `/usr/bin/python3` has neither. Anything that doesn't need xgboost can be tested
  locally: parsers, `refine_regions.py`, `modelling/src/events.py`.
- **Shell is zsh.** Never use `path` as a loop variable: it is tied to `$PATH`, and assigning
  it wipes command lookup. Unquoted globs such as `--include=*.py` error with "no matches
  found", so quote them.
- **Files copied from the external drive** arrive with macOS `._*` AppleDouble files and
  `rwx------` modes. `._*` is gitignored, but restore git's file modes after any bulk copy.

## Checks you can run locally

```bash
PY=~/miniforge3/envs/cnv_consensus_ml/bin/python
M=modules/xgboost_consensus_model

# row-to-event matching self-test (hand-checked expectations)
PYTHONDONTWRITEBYTECODE=1 $PY -c "import sys; sys.path.insert(0,'$M/modelling/src'); import events; events._self_test()"

# every Python file parses
$PY -c "import ast,pathlib; [ast.parse(p.read_text(),str(p)) for p in pathlib.Path('$M').rglob('*.py')]"

# workflow descriptor is valid JSON
python3 -I -c "import json; json.load(open('workflow.config.json'))"
```

On the cluster: `nextflow run main.nf -c local.config -stub -resume` exercises the full
channel wiring. Every process the workflow uses defines a `stub:` block.

## Invariants — breaking these fails silently

1. **Train/serve feature parity.** `FEATURE_COLUMNS` in `score_samples_xgb.py` must match
   `modelling/src/train_no_freec_cnmops.py` and the model's stored `feature_names`, in the
   same order. Change one and you must change all three, then retrain. Scoring with a
   mismatched order gives wrong probabilities, not an error.
2. **All 8 model callers stay enabled:** `cnvkit freec exomeDepth panelcnMOPS cnMOPS XHMM
   conifer gatk`. `freec` and `cnMOPS` are not direct features, but they feed
   `agreement_count` and candidate merging. Removing them from the scorer, or disabling them,
   shifts the inputs away from what the model was trained on. The `xgb_cutoff = 0.7499`
   operating point (precision ≈0.90) then stops meaning what the report says.
3. **The model never picks coordinates.** `_merge_candidates()` makes union coordinates, the
   model scores them, and `refine_regions.py` narrows the survivors to supported panel exons
   *before* `CLASSIFY_AND_ANNOTATE`. Don't move refinement after annotation. Annotating
   union spans attributes genes and classifications to places no caller called.
4. **Output schema.** Both combiners must emit `_smoothed_variants.tsv/.bed` in the same
   13-column `MERGE_VARIANT_CALLS` schema. Annotation and `FINAL_FORMATTING_TABLES` depend on
   it, and the reshape is the `awk` block at the end of the module's `main.nf`.
5. **Errors must fail the task.** `score_samples_xgb.py` writes a header-only table itself
   when a sample has no calls. Never reintroduce a `|| { write empty table }` fallback in
   `main.nf`: it turns real errors (missing xgboost, bad model) into "sample has no CNVs".
6. **The shipped model needs `xgboost >= 3.4`.** It was saved by 3.4.1, and
   `modelling/requirements.txt` pins the training environment. Keep `env.yaml` compatible.
7. **Sample names are full BAM basenames.** BAMs resolve to `mapped/<sample_name>.bam`, and
   `meta.sample_name` names every output. Ground truth must use the same names.

## Config: two sources of truth

- `nextflow.config` holds every param.
- `workflow.config.json` is the descriptor for the institutional workflow manager. Its
  `gui_params` defaults **override** `nextflow.config` for runs launched from that manager.
  When you add, rename or change the default of a user-facing param, update both, using the
  exact `params.<name>` the code reads. A name mismatch silently does nothing: the GUI
  toggle was once `use_control_freec` while the code read `use_controlfreec`.
- There is no float type in `gui_params`, so float params such as `xgb_cutoff` are not
  exposed there.
- `assembly` is functional, not informational:
  - ClassifyCNV uses hg19 for `GRCh37`/`hg19` and hg38 for anything else.
  - GATK uses `chrX` for `hg19`/`hg38` and `X` for anything else.
  - It must match the FASTA's build and chromosome naming.
- `min_number_callers` and `size_threshold` only affect the rule-based path.
  `smooth_gap_bp` is used by both paths.

## Known traps in the Nextflow code

- **Channel shape.** `Utils.mixAndCollectVarcalls()` groups by the first tuple element, so
  everything mixed into `ch_all_varcalls` must be `[meta, file(s)]`. ECOLE emits a bare
  cohort-level file list, so it is deliberately **not** mixed in. Doing so creates a fake
  sample, and `meta.sample_name` fails downstream. jabCoNtool's single cohort file is
  appended to every sample on purpose, and the model ignores it.
- **Process signature.** `XGBOOST_CONSENSUS_SCORE` takes 3 inputs: sample tuple, model JSON
  and panel BED. It has two call sites, `workflows/PANEL.nf` and
  `modelling/xgboost_pipeline/workflows/XGBOOST_SCORING_WES.nf`. Change both together.
- **Empty panel channel.** An unset `organism_dna_panel` gives an empty channel, and then the
  scorer never runs. Nothing reports an error.
- **Classic parser only.** The `.nf` files declare channels at script top level. That needs
  Nextflow's classic parser; the strict `v2` parser rejects it. Language-server errors like
  "Statements cannot be mixed with script declarations" come from this pre-existing layout,
  not from your edit.
- **`cache false`.** `XGBOOST_CONSENSUS_SCORE` and `FINAL_FORMATTING_TABLES` re-run on every
  `-resume`.
- **Decimal formatting.** Filenames built from floats in Groovy must use
  `String.format(java.util.Locale.ROOT, ...)` to match Python's `f"{x:.4f}"`. The cluster
  locale may be `cs_CZ`, which would write a decimal comma.

## Data protection

This is clinical genomics code. **Never commit patient data**, and never put real sample
identifiers into tracked files beyond what is already in `nextflow.config`.

- `modelling/xgboost_pipeline/test_prague.config` and `modelling/real_data_eval/results/`
  are gitignored because they hold real patient data. Keep them ignored. If you move them,
  move their `.gitignore` entries too.
- Real-sample test runs go in an untracked `-c` override file, never in a tracked config.
- Don't send repo data (calls, sample lists, BAM paths) to external services.

## Git conventions

- Branches: work on `dev_nextflow_model`. PRs target `plusWES`, the main branch.
- Commit messages use gitmoji plus a conventional-commit subject:
  `<emoji> <type>(<scope>): <imperative summary>`. For example:
  - `🐛 fix(events): credit a row with every ground-truth event it overlaps`
  - `✨ feat(xgboost-pipeline): add standalone XGBoost consensus-scoring pipeline`

  Types and emoji in use: ✨ feat, 🐛 fix, 📚 docs, ♻️ refactor, ✅ test, 🔨 chore, 📦 build,
  🔒 security.
- Commit only when asked. Keep one logical change per commit.
- `CLAUDE.md`, `IMPLEMENTATION_PLAN.md`, session-handoff and manuscript files are
  gitignored and stay local.

## Where to look

| Task | Start at |
|---|---|
| Add or toggle a caller | `workflows/PANEL.nf`, `subworkflows/<caller>/main.nf`, `nextflow.config` caller toggles, `workflow.config.json` |
| Change how calls are parsed for the model | `CALLER_FILE_PATTERNS` and the `parse_*` functions in `cnv_consensus_model.py` |
| Change scoring or the cutoff | `score_samples_xgb.py`, `xgb_cutoff` in `nextflow.config`, the module's `main.nf` |
| Change region refinement | `refine_regions.py` (`MAX_SUPPORT_GENES`, `smooth()`) |
| Retrain the model | `modelling/` steps in README "Training the model"; copy the new JSON to the module root |
| Cohort tables / HTML report | `bin/combine_final_tables.R`, `modules/final_formatting_tables/` |
| Annotation | `modules/classify_and_annotate_CNVs/` |
