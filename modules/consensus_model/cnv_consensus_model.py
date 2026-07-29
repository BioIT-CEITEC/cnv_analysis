"""
CNV Caller Consensus Model  (simulated-variant benchmark version)
=================================================================
Adapted from ~/Documents/cnv_analysis/training_data/
              structural_varcalls_3_callers_without_freec_and_cnvkit/cnv_consensus_model.py

What is different from the reference model
------------------------------------------
1. Ground truth are the SIMULATED variants (6 per sample) listed in
   selected_regions_simulated/<SAMPLE>_selected_regions.tsv.  They live in the
   *simulated* BAMs, so caller calls are read from `structural_varcalls_simulated_BR`.

2. "Baseline cleans the false positives" (user decision).
   Each simulated sample has a matching baseline sample (same name) in
   `structural_varcalls_baseline_BR` that was called WITHOUT the injected
   variants.  A simulated call that does not overlap any ground-truth variant
   is only counted as a false positive (negative training example) when the
   SAME caller did NOT already make that call in the baseline sample — i.e. it
   is a genuinely new call caused by the injection, not a pre-existing CNV.

3. Per-tool weights are constrained to the closed interval [0, weight_cap]
   (default weight_cap = 1.0) and can never be negative (user decision).
   The model is therefore NOT the reference StandardScaler + LogisticRegression
   pipeline (whose coefficients are signed and unbounded) but a small
   box-constrained logistic regression, `BoundedConsensusModel`, with ONE
   weight per caller and binary "did this caller detect the region" features.
   consensus_score(region) = sigmoid( intercept + Σ_c  weight_c * hit_c(region) )

4. Callers are restricted to the same set the reference model used
   (exomeDepth, panelcnMOPS, XHMM, conifer, gatk) via DEFAULT_CALLERS, but any
   subset present in the data can be selected on the command line.

The caller parsers, ground-truth loader and overlap logic are copied verbatim
from the reference so behaviour stays comparable.
"""

import json
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Every caller the parsers understand.
ALL_CALLERS = ["cnvkit", "freec", "exomeDepth", "panelcnMOPS", "cnMOPS",
               "XHMM", "conifer", "gatk"]

# The five callers the reference model was trained on ("los mismos callers").
DEFAULT_CALLERS = ["exomeDepth", "panelcnMOPS", "XHMM", "conifer", "gatk"]

# Minimum fraction of the SMALLER interval that must be covered for a match.
# (Kept identical to the reference: a large caller segment containing a small
#  GT exon still counts as a hit.)
MIN_OVERLAP_FRAC = 0.5

# When building merged candidates, any call whose size exceeds
# MAX_CLUSTER_SIZE_RATIO × the median call size in the cluster is an outlier
# (e.g. conifer reporting a gene-level deletion while exon-level callers report
# a single exon).  Outlier calls are counted in n_callers but their coordinates
# are excluded from the candidate boundaries so they cannot inflate a small
# exon-level candidate to a multi-megabase region.
MAX_CLUSTER_SIZE_RATIO = 50

# Sub-directories inside sv_dir that are not samples.
IGNORE_DIRS = {"all_samples", "selected_regions", "merged_variants",
               "selected_regions_simulated"}


# ---------------------------------------------------------------------------
# Utility helpers  (verbatim from reference)
# ---------------------------------------------------------------------------

def norm_chr(c):
    """Normalise chromosome name to 'chrN' form."""
    c = str(c).strip()
    if not c.startswith("chr"):
        c = "chr" + c
    return c


def norm_type(t):
    """Return 'DEL' or 'DUP' or None."""
    if pd.isna(t):
        return None
    t = str(t).lower()
    if any(k in t for k in ["del", "loss", "cn0", "cn1"]):
        return "DEL"
    if any(k in t for k in ["dup", "gain", "amp"]):
        return "DUP"
    m = re.match(r"cn(\d+)", t)
    if m:
        cn = int(m.group(1))
        if cn < 2:
            return "DEL"
        if cn > 2:
            return "DUP"
    return None


def gt_coverage(s1, e1, s2, e2):
    """Fraction of the SMALLER interval covered by its overlap with the larger."""
    overlap = max(0, min(e1, e2) - max(s1, s2))
    if overlap == 0:
        return 0.0
    min_len = min(e1 - s1, e2 - s2)
    if min_len <= 0:
        return 0.0
    return overlap / min_len


def overlaps(s1, e1, s2, e2, min_frac=0.0):
    """True if the two intervals overlap with min-length coverage >= min_frac."""
    return gt_coverage(s1, e1, s2, e2) >= min_frac


# ---------------------------------------------------------------------------
# Ground-truth loading  (Format B: *_selected_regions.tsv)
# ---------------------------------------------------------------------------

def merge_gt_variants(variants):
    """
    Deduplicate and merge exon-level GT entries into variant-level events.
    Two entries merge when on the same chr, same type, and start within 10 kb
    of the current end (consecutive panel exons of the same gene).
    """
    if not variants:
        return []
    GAP = 10_000
    seen, unique = set(), []
    for v in variants:
        key = (v["chr"], v["start"], v["end"], v["type"])
        if key not in seen:
            seen.add(key)
            unique.append(v)
    unique.sort(key=lambda x: (x["chr"], x["type"], x["start"]))
    merged = []
    cur = dict(unique[0])
    for v in unique[1:]:
        if (v["chr"] == cur["chr"] and v["type"] == cur["type"]
                and v["start"] <= cur["end"] + GAP):
            cur["end"] = max(cur["end"], v["end"])
        else:
            merged.append(cur)
            cur = dict(v)
    merged.append(cur)
    return merged


def load_ground_truth_dir(gt_dir):
    """
    Directory of per-sample TSVs named <SAMPLE>_selected_regions.tsv.
    Columns: sample, chrom, start, end, name, variant_type, type
    'type' holds DEL / DUP directly.
    Returns (gt: dict sample -> [ {chr,start,end,type} ], all_samples: set).
    """
    gt, all_samples = {}, set()
    for tsv in sorted(Path(gt_dir).glob("*_selected_regions.tsv")):
        if tsv.name.startswith(".") or tsv.stat().st_size == 0:
            continue
        df = pd.read_csv(tsv, sep="\t", encoding="utf-8", encoding_errors="replace")
        df.columns = df.columns.str.strip()
        for _, row in df.iterrows():
            sample = str(row["sample"]).strip()
            all_samples.add(sample)
            vtype = str(row["type"]).strip().upper()
            if vtype not in ("DEL", "DUP"):
                continue
            gt.setdefault(sample, []).append({
                "chr":   norm_chr(row["chrom"]),
                "start": int(row["start"]),
                "end":   int(row["end"]),
                "type":  vtype,
            })
    gt = {s: merge_gt_variants(v) for s, v in gt.items()}
    return gt, all_samples


def detect_callers(sv_dir):
    """
    Return callers (from ALL_CALLERS, in order) present in `sv_dir`.

    Two layouts are recognised, so detection matches what load_sample_calls()
    is actually able to read:
      * sv_dir/<sample>/<caller>/...        (caller-named sub-directory)
      * sv_dir/<sample>/<sample>_<caller>.* (flat, as staged by the pipeline)
    """
    found = set()
    for sample_dir in Path(sv_dir).iterdir():
        if not sample_dir.is_dir() or sample_dir.name in IGNORE_DIRS:
            continue
        for sub in sample_dir.iterdir():
            if sub.is_dir() and sub.name in ALL_CALLERS:
                found.add(sub.name)
        for caller in ALL_CALLERS:
            if caller in found:
                continue
            pattern = CALLER_FILE_PATTERNS.get(caller, f"*_{caller}.tsv")
            if next(sample_dir.rglob(pattern), None) is not None:
                found.add(caller)
            elif pattern.endswith(".vcf.gz") and next(sample_dir.rglob(pattern[:-3]), None) is not None:
                found.add(caller)
    return [c for c in ALL_CALLERS if c in found]


# ---------------------------------------------------------------------------
# Per-caller parsers  (verbatim from reference)
# Each returns a list of {chr, start, end, type, score}.
# ---------------------------------------------------------------------------

def parse_cnvkit(path):
    calls = []
    try:
        df = pd.read_csv(path, sep="\t")
        for _, r in df.iterrows():
            vtype = norm_type(r.get("SVTYPE", ""))
            if vtype is None:
                continue
            score = abs(float(r.get("FOLD_CHANGE_LOG", 0) or 0))
            calls.append({"chr": norm_chr(r["CHROM"]), "start": int(r["POS"]),
                          "end": int(r["END"]), "type": vtype, "score": score})
    except Exception:
        pass
    return calls


def parse_freec(path):
    calls = []
    try:
        df = pd.read_csv(path, sep="\t")
        for _, r in df.iterrows():
            vtype = norm_type(r.get("type", ""))
            if vtype is None:
                cn = r.get("CN", 2)
                try:
                    cn = float(cn)
                    vtype = "DEL" if cn < 2 else ("DUP" if cn > 2 else None)
                except Exception:
                    pass
            if vtype is None:
                continue
            cn = r.get("CN", 2)
            try:
                score = abs(float(cn) - 2.0)
            except Exception:
                score = 1.0
            calls.append({"chr": norm_chr(r["chromosome"]), "start": int(r["start"]),
                          "end": int(r["end"]), "type": vtype, "score": score})
    except Exception:
        pass
    return calls


def parse_exomeDepth(path):
    calls = []
    try:
        df = pd.read_csv(path, sep="\t")
        for _, r in df.iterrows():
            vtype = norm_type(r.get("type", ""))
            if vtype is None:
                continue
            bf = r.get("BF", 0)
            try:
                score = float(bf) if not pd.isna(bf) else 0.0
            except Exception:
                score = 0.0
            calls.append({"chr": norm_chr(r["chromosome"]), "start": int(r["start"]),
                          "end": int(r["end"]), "type": vtype, "score": score,
                          "cn": None})  # ExomeDepth does not report CN
    except Exception:
        pass
    return calls


def parse_panelcnMOPS(path):
    calls = []
    try:
        df = pd.read_csv(path, sep="\t")
        for _, r in df.iterrows():
            cn_raw = r.get("CN", "")
            vtype = norm_type(cn_raw)
            if vtype is None:
                continue
            try:
                rc_norm = float(r.get("RC.norm", 0) or 0)
                med_norm = float(r.get("medRC.norm", 1) or 1)
                score = abs(rc_norm / max(med_norm, 1) - 1)
            except Exception:
                score = 1.0
            # CN column is "CN0", "CN1", "CN3" etc. — extract the integer.
            cn = None
            try:
                cn_str = str(cn_raw).strip()
                if cn_str.upper().startswith("CN"):
                    cn = int(cn_str[2:])
            except Exception:
                pass
            calls.append({"chr": norm_chr(r["chromosome"]), "start": int(r["start"]),
                          "end": int(r["end"]), "type": vtype, "score": score, "cn": cn})
    except Exception:
        pass
    return calls


def parse_cnMOPS(path):
    calls = []
    try:
        df = pd.read_csv(path, sep="\t")
        if df.empty:
            return calls
        for _, r in df.iterrows():
            vtype = norm_type(r.get("CN", ""))
            if vtype is None:
                continue
            score = abs(float(r.get("mean", 0) or 0))
            calls.append({"chr": norm_chr(r["seqnames"]), "start": int(r["start"]),
                          "end": int(r["end"]), "type": vtype, "score": score})
    except Exception:
        pass
    return calls


def parse_XHMM(path):
    calls = []
    try:
        df = pd.read_csv(path, sep="\t")
        for _, r in df.iterrows():
            vtype = norm_type(r.get("type", ""))
            cn_raw = r.get("CN", 2)
            if vtype is None:
                try:
                    cn_raw = float(cn_raw)
                    vtype = "DEL" if cn_raw < 2 else ("DUP" if cn_raw > 2 else None)
                except Exception:
                    pass
            if vtype is None:
                continue
            cn_raw = r.get("CN", 2)
            cn = None
            try:
                cn = int(float(cn_raw))
            except Exception:
                pass
            try:
                score = abs(float(cn_raw) - 2.0)
            except Exception:
                score = 1.0
            calls.append({"chr": norm_chr(r["chromosome"]), "start": int(r["start"]),
                          "end": int(r["end"]), "type": vtype, "score": score, "cn": cn})
    except Exception:
        pass
    return calls


def parse_conifer(path):
    calls = []
    try:
        df = pd.read_csv(path, sep="\t")
        for _, r in df.iterrows():
            vtype = norm_type(r.get("type", "")) or norm_type(r.get("CN", ""))
            if vtype is None:
                continue
            try:
                rc_norm = float(r.get("RC.norm", 0) or 0)
                med_norm = float(r.get("medRC.norm", 1) or 1)
                score = abs(rc_norm / max(med_norm, 1) - 1)
            except Exception:
                score = 1.0
            calls.append({"chr": norm_chr(r["chromosome"]), "start": int(r["start"]),
                          "end": int(r["end"]), "type": vtype, "score": score,
                          "cn": None})  # conifer CN column is always empty
    except Exception:
        pass
    return calls


def parse_gatk(path):
    """GATK gCNV VCF (gzipped/BGZF or plain). ALT <DEL>/<DUP>, END in INFO, QS score."""
    import gzip
    import io as _io
    import subprocess
    import sys
    calls = []
    path = Path(path)
    if not path.exists():
        return calls
    if str(path).endswith(".vcf.gz") and path.stat().st_size == 4096:
        sib = Path(str(path)[:-3])
        if sib.exists():
            return parse_gatk(sib)

    raw_content = None
    for opener_fn in (lambda p: gzip.open(str(p), "rb").read(),
                      lambda p: open(str(p), "rb").read()):
        try:
            data = opener_fn(path)
            if data:
                raw_content = data
                break
        except Exception:
            pass
    if raw_content is None:
        for cmd in (["zcat", str(path)], ["gunzip", "-c", str(path)]):
            try:
                result = subprocess.run(cmd, stdout=subprocess.PIPE,
                                        stderr=subprocess.DEVNULL)
                if result.stdout:
                    raw_content = result.stdout
                    break
            except FileNotFoundError:
                pass
    if raw_content is None:
        print(f"WARNING: parse_gatk could not read {path}", file=sys.stderr)
        return calls

    header = None
    for raw in _io.BytesIO(raw_content):
        line = raw.decode("utf-8", errors="replace") if isinstance(raw, bytes) else raw
        if line.startswith("##"):
            continue
        if line.startswith("#"):
            header = line.lstrip("#").strip().split("\t")
            continue
        if header is None:
            continue
        row = dict(zip(header, line.strip().split("\t")))
        vtype = norm_type(row.get("ALT", "."))
        if vtype is None:
            continue
        chrom = norm_chr(row.get("CHROM", ""))
        pos = int(row.get("POS", 0))
        end = pos
        for token in row.get("INFO", "").split(";"):
            if token.startswith("END="):
                end = int(token.split("=")[1])
                break
        fmt_keys = row.get("FORMAT", "").split(":")
        sample_col = [v for k, v in row.items()
                      if k not in ("CHROM", "POS", "ID", "REF", "ALT",
                                   "QUAL", "FILTER", "INFO", "FORMAT")]
        score = 0.0
        cn = None
        if sample_col:
            fmt = dict(zip(fmt_keys, sample_col[0].split(":")))
            try:
                score = float(fmt.get("QS", fmt.get("QA", 0)) or 0)
            except (ValueError, TypeError):
                score = 0.0
            try:
                cn = int(fmt["CN"])
            except (KeyError, ValueError, TypeError):
                cn = None
        calls.append({"chr": chrom, "start": pos, "end": end,
                      "type": vtype, "score": score, "cn": cn})
    return calls


CALLER_PARSERS = {
    "cnvkit": parse_cnvkit, "freec": parse_freec, "exomeDepth": parse_exomeDepth,
    "panelcnMOPS": parse_panelcnMOPS, "cnMOPS": parse_cnMOPS, "XHMM": parse_XHMM,
    "conifer": parse_conifer, "gatk": parse_gatk,
}

CALLER_FILE_PATTERNS = {
    "cnvkit": "*_cnvkit.tsv", "freec": "*_freec.tsv", "exomeDepth": "*_ExomeDepth.tsv",
    "panelcnMOPS": "*_panelcnMOPS.tsv", "cnMOPS": "*_cnMOPS.tsv", "XHMM": "*_xhmm.tsv",
    "conifer": "*_conifer.tsv", "gatk": "*_gatk.vcf.gz",
}


def load_sample_calls(sample_dir, callers):
    """Return {caller: [call dicts]} for one sample directory."""
    calls_by_caller = {}
    for caller in callers:
        pattern = CALLER_FILE_PATTERNS.get(caller, f"*_{caller}.tsv")
        matches = list(Path(sample_dir).rglob(pattern))
        if not matches and pattern.endswith(".vcf.gz"):
            matches = list(Path(sample_dir).rglob(pattern[:-3]))
        if not matches:
            calls_by_caller[caller] = []
            continue
        best = next((m for m in matches if m.parent.name == caller), matches[0])
        calls_by_caller[caller] = CALLER_PARSERS[caller](best)
    return calls_by_caller


# ---------------------------------------------------------------------------
# Matching helpers
# ---------------------------------------------------------------------------

def _best_hit(target, calls):
    """(hit_binary, best_score) for `calls` vs `target` region (same chr+type, overlap)."""
    best_score, hit = 0.0, 0
    for c in calls:
        if c["chr"] != target["chr"]:
            continue
        if c.get("type") and target.get("type") and c["type"] != target["type"]:
            continue
        if overlaps(c["start"], c["end"], target["start"], target["end"], MIN_OVERLAP_FRAC):
            hit = 1
            best_score = max(best_score, c.get("score", 0.0))
    return hit, best_score


# ---------------------------------------------------------------------------
# Feature matrix  (reference strategy: [hit, score] per caller, simulated only)
# ---------------------------------------------------------------------------

def build_feature_matrix(sim_dir, gt, samples, callers):
    """
    Build the feature matrix exactly as in the reference cnv_analysis model:
    two columns per caller, [hit, score].

      hit   = 1 if the caller detected the candidate region (>= MIN_OVERLAP_FRAC)
      score = the caller's confidence for its best matching call (0 if no hit)

    Positives : one row per ground-truth variant (features from SIMULATED calls).
    Negatives : each SIMULATED call that does not overlap any GT variant
                (raw / genome-wide false positives).  No baseline is used.

    Returns (X, y, feature_names, meta).
    """
    rows, labels, meta = [], [], []

    for sample in sorted(samples):
        sim_sd = Path(sim_dir) / sample
        if not sim_sd.is_dir():
            continue
        sample_gt = gt.get(sample, [])
        sim_calls = load_sample_calls(sim_sd, callers)

        # ---- Positives ----
        for g in sample_gt:
            row = []
            for c in callers:
                hit, score = _best_hit(g, sim_calls.get(c, []))
                row.extend([hit, score])
            rows.append(row)
            labels.append(1)
            meta.append((sample, g["chr"], g["start"], g["end"], g["type"], "TP"))

        # ---- Negatives ----
        for caller in callers:
            for call in sim_calls.get(caller, []):
                is_tp = any(call["chr"] == g["chr"] and call["type"] == g["type"]
                            and overlaps(call["start"], call["end"],
                                         g["start"], g["end"], MIN_OVERLAP_FRAC)
                            for g in sample_gt)
                if is_tp:
                    continue
                row = []
                for c2 in callers:
                    if c2 == caller:
                        row.extend([1, call.get("score", 0.0)])
                    else:
                        hit2, score2 = _best_hit(call, sim_calls.get(c2, []))
                        row.extend([hit2, score2])
                rows.append(row)
                labels.append(0)
                meta.append((sample, call["chr"], call["start"], call["end"],
                             call["type"], "FP"))

    feature_names = []
    for c in callers:
        feature_names.extend([f"{c}_hit", f"{c}_score"])
    X = np.array(rows, dtype=float) if rows else np.zeros((0, 2 * len(callers)))
    y = np.array(labels, dtype=int)
    return X, y, feature_names, meta


# ---------------------------------------------------------------------------
# Consensus model  (reference strategy + non-negative capped coefficients)
# ---------------------------------------------------------------------------

class BoundedConsensusModel:
    """
    The reference cnv_analysis consensus model — StandardScaler followed by a
    logistic regression on [hit, score] features — with ONE change: every
    coefficient is constrained to be NON-NEGATIVE and capped at `weight_cap`
    (default 1.0).  So a caller's presence/confidence can only ever *increase*
    the consensus probability, never decrease it, and no single caller can
    dominate beyond the cap.

        p = sigmoid( intercept + Σ_f coef_f * z_f ),   z = (x - mean) / scale
        coef_f ∈ [0, weight_cap]     (intercept is free)

    Coefficients live in the standardised feature space, exactly like the
    reference `hit_coef` / `score_coef`.  Fitted by L-BFGS-B on the
    class-weighted logistic loss with a small L2 penalty.  Picklable and
    self-contained (stores the scaler + the caller list for scoring).
    """

    def __init__(self, feature_names, callers=None, weight_cap=1.0, l2=0.01,
                 class_weight="balanced", cap_mode="per_caller"):
        self.feature_names = list(feature_names)
        self.callers = list(callers) if callers is not None else None
        self.weight_cap = float(weight_cap)
        self.l2 = float(l2)
        self.class_weight = class_weight
        # cap_mode controls what "no weight exceeds weight_cap" means:
        #   "per_caller" -> each caller's combined weight (hit+score) <= weight_cap
        #                   (guarantees every reported number <= weight_cap)
        #   "per_coef"   -> each individual coefficient <= weight_cap (combined
        #                   can reach 2*weight_cap)
        self.cap_mode = cap_mode
        self.coef_ = None          # coefficients in standardised space
        self.intercept_ = 0.0
        self.mean_ = None          # StandardScaler params
        self.scale_ = None

    # -- internals -------------------------------------------------------
    @staticmethod
    def _sigmoid(z):
        return np.where(z >= 0, 1.0 / (1.0 + np.exp(-z)),
                        np.exp(z) / (1.0 + np.exp(z)))

    def _sample_weights(self, y):
        sw = np.ones_like(y, dtype=float)
        if self.class_weight == "balanced":
            n = len(y)
            n_pos = max(int(y.sum()), 1)
            n_neg = max(int((y == 0).sum()), 1)
            sw[y == 1] = n / (2.0 * n_pos)
            sw[y == 0] = n / (2.0 * n_neg)
        return sw

    def _scale(self, X):
        return (np.asarray(X, dtype=float) - self.mean_) / self.scale_

    # -- API -------------------------------------------------------------
    def fit(self, X, y):
        from scipy.optimize import minimize
        X = np.asarray(X, dtype=float)
        y = np.asarray(y, dtype=int)

        # StandardScaler (guard constant columns like the reference does).
        self.mean_ = X.mean(axis=0)
        scale = X.std(axis=0)
        scale[scale == 0] = 1.0
        self.scale_ = scale
        Xs = self._scale(X)

        k = Xs.shape[1]
        sw = self._sample_weights(y)
        sw_sum = sw.sum()

        def objective(params):
            w = params[:k]
            b = params[k]
            z = Xs @ w + b
            loss = np.sum(sw * (np.logaddexp(0, z) - y * z)) / sw_sum
            loss += self.l2 * np.dot(w, w)
            resid = sw * (self._sigmoid(z) - y)
            grad_w = Xs.T @ resid / sw_sum + 2.0 * self.l2 * w
            grad_b = resid.sum() / sw_sum
            return loss, np.append(grad_w, grad_b)

        bounds = [(0.0, self.weight_cap)] * k + [(None, None)]

        # Per-caller cap: constrain the SUM of a caller's coefficients
        # (hit + score) to <= weight_cap, so no reported weight ever exceeds it.
        groups = []
        if self.cap_mode == "per_caller" and self.callers:
            for c in self.callers:
                idx = [i for i, f in enumerate(self.feature_names)
                       if f in (f"{c}_hit", f"{c}_score")]
                if idx:
                    groups.append(idx)

        if groups:
            constraints = []
            for idx in groups:
                jac = np.zeros(k + 1)
                jac[idx] = -1.0
                constraints.append({
                    "type": "ineq",
                    "fun": (lambda p, idx=idx: self.weight_cap - float(np.sum(p[idx]))),
                    "jac": (lambda p, jac=jac: jac),
                })
            x0 = np.append(np.full(k, min(0.3, self.weight_cap)), 0.0)
            res = minimize(objective, x0, jac=True, method="SLSQP",
                           bounds=bounds, constraints=constraints,
                           options={"maxiter": 1000, "ftol": 1e-8})
        else:
            x0 = np.append(np.full(k, min(0.5, self.weight_cap)), 0.0)
            res = minimize(objective, x0, jac=True, method="L-BFGS-B", bounds=bounds)

        self.coef_ = res.x[:k]
        self.intercept_ = float(res.x[k])
        return self

    def decision_function(self, X):
        return self._scale(X) @ self.coef_ + self.intercept_

    def predict_proba(self, X):
        p = self._sigmoid(self.decision_function(X))
        return np.column_stack([1.0 - p, p])

    def feature_weight_dict(self):
        """{feature_name: coef}."""
        return {f: round(float(w), 4)
                for f, w in zip(self.feature_names, self.coef_)}

    def caller_weight_dict(self):
        """{caller: {hit_coef, score_coef, combined}} — same shape as reference."""
        fw = dict(zip(self.feature_names, self.coef_))
        callers = self.callers or [f[:-4] for f in self.feature_names
                                   if f.endswith("_hit")]
        out = {}
        for c in callers:
            hit = round(float(fw.get(f"{c}_hit", 0.0)), 4)
            score = round(float(fw.get(f"{c}_score", 0.0)), 4)
            out[c] = {"hit_coef": hit, "score_coef": score,
                      "combined": round(hit + score, 4)}
        return out


# ---------------------------------------------------------------------------
# Scoring new samples
# ---------------------------------------------------------------------------

def _merge_candidates(df):
    """Merge overlapping calls (same chr+type) into representative regions.

    After the greedy union merge, calls whose size exceeds MAX_CLUSTER_SIZE_RATIO
    times the median cluster call size are treated as coordinate outliers: they
    are counted in n_callers / callers_supporting but their boundaries are
    excluded from the final candidate coordinates.  This prevents a single
    gene-level caller (e.g. conifer) from inflating a candidate that is
    otherwise supported by exon-level callers at a much smaller scale.
    """
    merged_rows = []
    for (chrom, vtype), group in df.groupby(["chr", "type"]):
        group = group.sort_values("start").reset_index(drop=True)
        used = [False] * len(group)
        for i in range(len(group)):
            if used[i]:
                continue
            caller_i = group.at[i, "caller"]
            cluster = {caller_i}
            cn_by_caller = {caller_i: group.at[i, "cn"] if "cn" in group.columns else None}
            cstart, cend = group.at[i, "start"], group.at[i, "end"]
            intervals = [(group.at[i, "start"], group.at[i, "end"])]
            for j in range(i + 1, len(group)):
                if used[j]:
                    continue
                if overlaps(cstart, cend, group.at[j, "start"],
                            group.at[j, "end"], MIN_OVERLAP_FRAC):
                    caller_j = group.at[j, "caller"]
                    cluster.add(caller_j)
                    if caller_j not in cn_by_caller:
                        cn_by_caller[caller_j] = group.at[j, "cn"] if "cn" in group.columns else None
                    cstart = min(cstart, group.at[j, "start"])
                    cend = max(cend, group.at[j, "end"])
                    intervals.append((group.at[j, "start"], group.at[j, "end"]))
                    used[j] = True
            used[i] = True

            # Clip outlier boundaries when cluster has more than one call.
            if len(intervals) > 1:
                sizes = sorted(e - s for s, e in intervals)
                median_size = sizes[len(sizes) // 2]
                threshold = median_size * MAX_CLUSTER_SIZE_RATIO
                small = [(s, e) for s, e in intervals if e - s <= threshold]
                if small and len(small) < len(intervals):
                    cstart = min(s for s, e in small)
                    cend = max(e for s, e in small)

            cn_label = ";".join(
                "{}:{}".format(c, cn_by_caller.get(c) if cn_by_caller.get(c) is not None else "NA")
                for c in sorted(cluster)
            )
            merged_rows.append({"chr": chrom, "start": cstart, "end": cend,
                                "type": vtype, "n_callers": len(cluster),
                                "callers_supporting": ",".join(sorted(cluster)),
                                "cn_label": cn_label})
    return pd.DataFrame(merged_rows)


def score_calls(calls_by_caller, model, callers=None):
    """
    Score one sample. Returns a DataFrame with columns:
      chr, start, end, type, n_callers, callers_supporting, consensus_score
    sorted by descending consensus_score.

    Features (binary hits, or hit+quality) are rebuilt exactly as at training
    time using the metadata stored on the model.
    """
    callers = callers or model.callers
    candidates = []
    for caller, calls in calls_by_caller.items():
        for call in calls:
            candidates.append({"chr": call["chr"], "start": call["start"],
                               "end": call["end"], "type": call["type"],
                               "caller": caller, "score": call.get("score", 0.0),
                               "cn": call.get("cn")})
    if not candidates:
        return pd.DataFrame()

    merged = _merge_candidates(pd.DataFrame(candidates))
    rows = []
    for _, reg in merged.iterrows():
        target = {"chr": reg["chr"], "start": reg["start"],
                  "end": reg["end"], "type": reg["type"]}
        row = []
        for c in callers:
            hit, score = _best_hit(target, calls_by_caller.get(c, []))
            row.extend([hit, score])
        rows.append(row)
    X_new = np.array(rows, dtype=float)
    merged = merged.copy()
    merged["consensus_score"] = model.predict_proba(X_new)[:, 1]
    # Keep consensus_score immediately after callers_supporting (col 8 in TSV),
    # with cn_label as the last column so the awk reformat in main.nf stays stable.
    col_order = ["chr", "start", "end", "type", "n_callers", "callers_supporting",
                 "consensus_score", "cn_label"]
    for col in col_order:
        if col not in merged.columns:
            merged[col] = "NA"
    return merged[col_order].sort_values("consensus_score", ascending=False).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Persistence
# ---------------------------------------------------------------------------

def load_pipeline(pkl_path):
    import pickle
    with open(pkl_path, "rb") as fh:
        return pickle.load(fh)
