"""Shared paths, covariate sets and estimators.

Every random step takes an explicit numpy Generator so that the scripts
reproduce the estimates reported in the manuscript exactly.
"""
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier as HGB
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold, cross_val_predict

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"          # outputs of the SQL queries (not distributed)
DERIVED = ROOT / "data" / "derived"  # analytic cohorts (not distributed)
RESULTS = ROOT / "results"           # aggregate results only
for p in (DERIVED, RESULTS):
    p.mkdir(parents=True, exist_ok=True)

THRESHOLD = 0.05  # NEE at discontinuation (ug/kg/min) separating abrupt from tapered stops
MAX_STOP_DOSE = 0.2
B = 500           # bootstrap resamples

# Covariates for the propensity score (MIMIC-IV); 'cu' = first care unit (categorical)
CONF_MIMIC = ["age", "male", "sofa_d1", "sepsis3", "icu_h_at_stop", "ep_dur_h", "max_d",
              "vaso_any", "vaso_late", "other_late", "inotrope_any", "map_m2", "mapmin_m2",
              "map_24", "frac_map65_24", "hr_m2", "hr_24", "lac_last", "lac_max", "uo24", "imv",
              "cr_last", "ph_last", "cu", "dnr_by_stop", "propofol_on", "fentanyl_on",
              "dexmed_on", "hydrocort"]
# eICU-CRD; 'ut' = unit type (categorical)
CONF_EICU = ["age", "male", "apache", "apache_pred", "icu_h_at_stop", "ep_dur_h", "max_d",
             "vaso_any", "vaso_late", "n_agents_late", "map_m2", "map_24", "frac_map65_24",
             "hr_m2", "hr_24", "lac_last", "lac_max", "cr_last", "ph_last", "imv",
             "dnr_by_stop", "ut"]


def load_mimic():
    return pd.read_csv(DERIVED / "mimic_analytic.csv", float_precision="round_trip")


def load_eicu():
    return pd.read_csv(DERIVED / "eicu_analytic.csv", float_precision="round_trip")


def cboot(groups, rng, n_boot=B):
    """Cluster bootstrap: resample clusters (patients or hospitals) with replacement."""
    ug = np.unique(groups)
    pos = {g: np.where(groups == g)[0] for g in ug}
    for _ in range(n_boot):
        s = rng.choice(ug, len(ug))
        yield np.concatenate([pos[g] for g in s])


def wrisk(t, w, y, i):
    """Weighted risk in exposed and unexposed, restricted to index i."""
    return (np.sum(w[i] * t[i] * y[i]) / np.sum(w[i] * t[i]),
            np.sum(w[i] * (1 - t[i]) * y[i]) / np.sum(w[i] * (1 - t[i])))


def est(t, w, y, groups, label, out, rng, n_boot=B):
    """Weighted risks, risk difference and risk ratio with cluster-bootstrap 95% CI."""
    y = np.asarray(y).astype(int)
    r1, r0 = wrisk(t, w, y, np.arange(len(t)))
    bs = np.array([wrisk(t, w, y, i) for i in cboot(groups, rng, n_boot)])
    rd = bs[:, 0] - bs[:, 1]
    rr = bs[:, 0] / np.maximum(bs[:, 1], 1e-9)
    res = dict(analysis=label, n=len(t), n_exp=int(t.sum()), r1=r1, r0=r0, rd=r1 - r0,
               rd_lo=np.percentile(rd, 2.5), rd_hi=np.percentile(rd, 97.5), rr=r1 / r0,
               rr_lo=np.percentile(rr, 2.5), rr_hi=np.percentile(rr, 97.5))
    out.append(res)
    print(f"{label:42s} n={len(t):>6} exp={int(t.sum()):>5}  {100*r1:5.1f} vs {100*r0:5.1f}  "
          f"RD {100*res['rd']:.1f} ({100*res['rd_lo']:.1f} to {100*res['rd_hi']:.1f})  "
          f"RR {res['rr']:.2f} ({res['rr_lo']:.2f}-{res['rr_hi']:.2f})")
    return res


def design(D, cols, cat):
    """Standardised design matrix: median-imputed covariates plus care-unit dummies."""
    X = D[cols].astype(float)
    X = X.fillna(X.median())
    X = pd.concat([X, pd.get_dummies(D[cat]).astype(float)], axis=1).values
    return (X - X.mean(0)) / (X.std(0) + 1e-9)


def aipw(t, X, y, ps):
    """Augmented inverse probability weighted risk difference."""
    o1 = LogisticRegression(max_iter=3000, C=0.5).fit(X[t == 1], y[t == 1]).predict_proba(X)[:, 1]
    o0 = LogisticRegression(max_iter=3000, C=0.5).fit(X[t == 0], y[t == 0]).predict_proba(X)[:, 1]
    return np.mean(o1 - o0 + t * (y - o1) / ps - (1 - t) * (y - o0) / (1 - ps))


def smdmax(D, cols, t, w):
    """Maximum absolute weighted standardised mean difference."""
    m = 0
    for c in cols:
        x = D[c].astype(float)
        x = x.fillna(x.median()).values
        s = np.sqrt((x[t == 1].var() + x[t == 0].var()) / 2)
        if s > 0:
            m = max(m, abs((np.average(x[t == 1], weights=w[t == 1])
                            - np.average(x[t == 0], weights=w[t == 0])) / s))
    return m


def match(t, ps, cal_sd=0.2, seed=1):
    """1:1 nearest-neighbour matching on the logit of the PS without replacement, caliper 0.2 SD."""
    r = np.random.default_rng(seed)
    lp = np.log(ps / (1 - ps))
    cal = cal_sd * lp.std()
    ct = np.where(t == 0)[0]
    o = np.argsort(lp[ct])
    cs, ci = lp[ct][o], ct[o]
    used = np.zeros(len(t), bool)
    pairs = []
    for i in r.permutation(np.where(t == 1)[0]):
        j = np.searchsorted(cs, lp[i])
        best = None
        for k in range(max(0, j - 30), min(len(cs), j + 30)):
            c = ci[k]
            if used[c]:
                continue
            dd = abs(cs[k] - lp[i])
            if dd < cal and (best is None or dd < best[0]):
                best = (dd, c)
        if best:
            used[best[1]] = True
            pairs.append((i, best[1]))
    return np.array(pairs)


def psfit(D, thr, groups, cols, catname):
    """Exposure indicator and cross-fitted gradient-boosted propensity score (clipped 0.01-0.99)."""
    t = (D.last_d > thr).astype(int).values
    model = HGB(max_iter=300, learning_rate=0.05,
                categorical_features=[cols.index(catname)], random_state=0)
    ps = cross_val_predict(model, D[cols], t, cv=GroupKFold(5), groups=groups,
                           method="predict_proba")[:, 1]
    return t, np.clip(ps, 0.01, 0.99)


def evalue(rr):
    return rr + np.sqrt(rr * (rr - 1))
