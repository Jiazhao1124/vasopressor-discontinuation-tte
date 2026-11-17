"""Primary, secondary, sensitivity and subgroup analyses (Table 2, Table S3, Figure S2).

Output: results/main_results.csv
"""
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from common import (CONF_MIMIC, CONF_EICU, RESULTS, B, load_mimic, load_eicu, est, cboot,
                    design, aipw, smdmax, match, psfit, evalue)

rng = np.random.default_rng(2026)
OUT = []

# ---------------------------------------------------------------- MIMIC-IV
d = load_mimic()
g = d.subject_id.values
num = [c for c in CONF_MIMIC if c != "cu"]
t, ps = psfit(d, 0.05, g, CONF_MIMIC, "cu")
w = np.where(t == 1, 1 - ps, ps)                      # overlap weights
X = design(d, num, "cu")
print("MIMIC-IV n", len(d), "patients", len(np.unique(g)),
      "max weighted SMD", round(smdmax(d, num, t, w), 3))

y24 = (d.restart_gap_h <= 24).astype(int).values
prim = est(t, w, y24, g, "M primary OW", OUT, rng)
est(t, w, (d.restart_gap_h > 2) & (d.restart_gap_h <= 24), g, "M restart >2h", OUT, rng)
est(t, w, d.d24nr == 1, g, "M death24 no restart", OUT, rng)
est(t, w, (d.restart_gap_h <= 24) | (d.death_h <= 24), g, "M composite", OUT, rng)
est(t, w, d.death28 == 1, g, "M death28", OUT, rng)
print("E-value: point", round(evalue(prim["rr"]), 2), "lower CI", round(evalue(prim["rr_lo"]), 2))

# IPTW, trimmed at PS 0.05-0.95
m = (ps > 0.05) & (ps < 0.95)
wi = np.where(t == 1, 1 / ps, 1 / (1 - ps))
est(t[m], wi[m], y24[m], g[m], "M IPTW trimmed", OUT, rng)

# The reported run bootstrapped an exploratory AIPW with the gradient-boosted PS here
# (unstable in eICU-CRD, not reported). Drawing the same 200 cluster resamples keeps the
# random stream, and therefore every later estimate, identical to the manuscript.
for _ in cboot(g, rng, 200):
    pass

# 1:1 propensity score matching, pair bootstrap
pr = match(t, ps)
dif = y24[pr[:, 0]] - y24[pr[:, 1]]
mb = [dif[rng.integers(0, len(dif), len(dif))].mean() for _ in range(B)]
OUT.append(dict(analysis="M PS-matched", n=2 * len(pr), n_exp=len(pr), rd=dif.mean(),
                rd_lo=np.percentile(mb, 2.5), rd_hi=np.percentile(mb, 97.5)))
print("M matched pairs", len(pr), round(100 * dif.mean(), 1), np.round(100 * np.percentile(mb, [2.5, 97.5]), 1))

# overlap weights from a logistic propensity model
pl = LogisticRegression(max_iter=3000, C=1.0).fit(X, t).predict_proba(X)[:, 1]
wl = np.where(t == 1, 1 - pl, pl)
print("M logistic-OW max weighted SMD", round(smdmax(d, num, t, wl), 3))
est(t, wl, y24, g, "M OW logistic PS", OUT, rng)

# first ICU stay per patient
fs = ~d.duplicated("subject_id", keep="first").values
est(t[fs], w[fs], y24[fs], g[fs], "M first stay per patient", OUT, rng)

# alternative thresholds for an abrupt stop
for thr in [0.02, 0.03, 0.08, 0.10]:
    tt, pp = psfit(d, thr, g, CONF_MIMIC, "cu")
    est(tt, np.where(tt == 1, 1 - pp, pp), y24, g, f"M threshold {thr:.2f}", OUT, rng)

cardiac = d.first_careunit.str.contains("Cardiac")
subs = [("Sepsis-3", d.sepsis3 == 1), ("No sepsis-3", d.sepsis3 == 0),
        ("Cardiovascular ICU", cardiac), ("Other ICUs", ~cardiac),
        ("Vasopressin at stop", d.vaso_late == 1), ("No vasopressin at stop", d.vaso_late == 0),
        ("Single agent, last 2 h", (d.other_late == 0) & (d.vaso_late == 0)),
        ("Episode >=24 h", d.ep_dur_h >= 24), ("Episode <24 h", d.ep_dur_h < 24),
        ("IMV at stop", d.imv == 1), ("No IMV at stop", d.imv == 0)]
for lab, mm in subs:
    mm = mm.values
    est(t[mm], w[mm], y24[mm], g[mm], "M sub: " + lab, OUT, rng, n_boot=300)

# ---------------------------------------------------------------- eICU-CRD
e = load_eicu()
num_e = [c for c in CONF_EICU if c != "ut"]
te, pe = psfit(e, 0.05, e.hospitalid.values, CONF_EICU, "ut")
we = np.where(te == 1, 1 - pe, pe)
ge, gh = e.patientunitstayid.values, e.hospitalid.values
Xe = design(e, num_e, "unittype")
ye = (e.restart_gap_h <= 24).astype(int).values
print("eICU-CRD max weighted SMD (HGB)", round(smdmax(e, num_e, te, we), 3))
est(te, we, ye, ge, "E primary OW", OUT, rng)
est(te, we, ye, gh, "E primary OW, hospital-cluster bootstrap", OUT, rng)
est(te, we, e.d24nr == 1, ge, "E death24 no restart", OUT, rng)
est(te, we, (e.restart_gap_h <= 24) | (e.hosp_death_h <= 24), ge, "E composite", OUT, rng)
est(te, we, (e.hosp_death == 1) & (e.hosp_death_h <= 672), ge, "E death28", OUT, rng)
pl = LogisticRegression(max_iter=3000).fit(Xe, te).predict_proba(Xe)[:, 1]
wl = np.where(te == 1, 1 - pl, pl)
print("eICU-CRD logistic-OW max weighted SMD", round(smdmax(e, num_e, te, wl), 3))
est(te, wl, ye, ge, "E OW logistic PS", OUT, rng)

# ---------------------------------------------------------------- AIPW (doubly robust)
# Logistic propensity score truncated at 0.05-0.95 and logistic outcome models, both refitted
# in each cluster-bootstrap resample.
def aipw_logistic(X, t, y, idx):
    p = np.clip(LogisticRegression(max_iter=3000).fit(X[idx], t[idx]).predict_proba(X[idx])[:, 1], 0.05, 0.95)
    return aipw(t[idx], X[idx], y[idx], p)


for lab, Xa, ta, ya, ga, nb in [("M AIPW", X, t, y24, g, 200), ("E AIPW", Xe, te, ye, ge, 300)]:
    r_ = np.random.default_rng(2026)
    a0 = aipw_logistic(Xa, ta, ya, np.arange(len(ta)))
    ab = [aipw_logistic(Xa, ta, ya, i) for i in cboot(ga, r_, nb)]
    OUT.append(dict(analysis=lab, n=len(ta), n_exp=int(ta.sum()), rd=a0,
                    rd_lo=np.percentile(ab, 2.5), rd_hi=np.percentile(ab, 97.5)))
    print(lab, round(100 * a0, 1), np.round(100 * np.percentile(ab, [2.5, 97.5]), 1))

pd.DataFrame(OUT).to_csv(RESULTS / "main_results.csv", index=False)
