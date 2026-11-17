"""Table S4: covariate balance before and after 1:1 propensity score matching.

Output: results/tableS4_psm_balance.csv, results/psm_matched_effects.csv
"""
import numpy as np
import pandas as pd
from common import CONF_MIMIC, CONF_EICU, RESULTS, load_mimic, load_eicu, psfit, match

d, e = load_mimic(), load_eicu()
t, ps = psfit(d, 0.05, d.subject_id.values, CONF_MIMIC, "cu"); pr = match(t, ps)
te, pe = psfit(e, 0.05, e.hospitalid.values, CONF_EICU, "ut"); pre = match(te, pe)

rng = np.random.default_rng(2026)
eff = []
for D, P, lab in [(d, pr, "MIMIC-IV"), (e, pre, "eICU-CRD")]:
    y = (D.restart_gap_h <= 24).astype(int).values
    dif = y[P[:, 0]] - y[P[:, 1]]
    mb = [dif[rng.integers(0, len(dif), len(dif))].mean() for _ in range(500)]
    eff.append(dict(database=lab, pairs=len(P), risk_abrupt=y[P[:, 0]].mean(), risk_tapered=y[P[:, 1]].mean(),
                    rd=dif.mean(), rd_lo=np.percentile(mb, 2.5), rd_hi=np.percentile(mb, 97.5)))
    print(eff[-1])
pd.DataFrame(eff).to_csv(RESULTS / "psm_matched_effects.csv", index=False)

CU_M = {"Cardiac Vascular Intensive Care Unit (CVICU)": "card", "Coronary Care Unit (CCU)": "card",
        "Medical Intensive Care Unit (MICU)": "med", "Medical/Surgical Intensive Care Unit (MICU/SICU)": "mix",
        "Surgical Intensive Care Unit (SICU)": "surg", "Trauma SICU (TSICU)": "surg"}
CU_E = {"Med-Surg ICU": "mix", "CCU-CTICU": "card", "Cardiac ICU": "card", "CTICU": "card",
        "CSICU": "card", "MICU": "med", "SICU": "surg"}
d["cu5"] = d.first_careunit.map(CU_M).fillna("oth"); e["cu5"] = e.unittype.map(CU_E).fillna("oth")
for D in (d, e):
    for k in ["card", "med", "surg", "mix", "oth"]:
        D["cu_" + k] = (D.cu5 == k).astype(int)
    D["frac65"] = D.frac_map65_24 * 100
d["apache_pred_pct"] = np.nan; e["apache_pred_pct"] = e.apache_pred * 100


def summ(x, kind, dec):
    x = x.dropna()
    if kind == "b": return f"{int(x.sum()):,} ({100*x.mean():.1f})"
    f = f"{{:,.{dec}f}}"; return f"{f.format(x.median())} ({f.format(x.quantile(.25))}–{f.format(x.quantile(.75))})"


def smd(x1, x0, kind):
    x1, x0 = x1.dropna().astype(float), x0.dropna().astype(float)
    if kind == "b":
        p1, p0 = x1.mean(), x0.mean(); s = np.sqrt((p1 * (1 - p1) + p0 * (1 - p0)) / 2)
    else:
        s = np.sqrt((x1.var() + x0.var()) / 2)
    return (x1.mean() - x0.mean()) / s if s > 0 else 0.0


V = [("age", "Age, years", "c", 1), ("male", "Male sex", "b", 0), ("sofa_d1", "SOFA score, ICU day 1", "c", 1),
     ("sepsis3", "Sepsis-3", "b", 0), ("apache", "APACHE IVa score", "c", 1),
     ("apache_pred_pct", "APACHE IVa predicted hospital mortality, %", "c", 1),
     ("cu_card", "Cardiac or cardiovascular ICU", "b", 0), ("cu_med", "Medical ICU", "b", 0),
     ("cu_surg", "Surgical or trauma ICU", "b", 0), ("cu_mix", "Mixed medical–surgical ICU", "b", 0),
     ("cu_oth", "Neurological or other ICU", "b", 0), ("ep_dur_h", "Episode duration, h", "c", 1),
     ("icu_h_at_stop", "Time from ICU admission to discontinuation, h", "c", 1), ("max_d", "Peak NEE, µg/kg/min", "c", 3),
     ("last_d", "NEE at discontinuation, µg/kg/min (defines exposure)", "c", 3),
     ("vaso_late", "Vasopressin at discontinuation", "b", 0), ("map_m2", "Mean arterial pressure, mmHg", "c", 1),
     ("hr_m2", "Heart rate, beats/min", "c", 1), ("frac65", "MAP readings <65 mmHg in preceding 24 h, %", "c", 1),
     ("lac_last", "Lactate, mmol/L", "c", 1), ("lac_max", "Peak lactate, mmol/L", "c", 1),
     ("cr_last", "Creatinine, mg/dL", "c", 1), ("ph_last", "pH", "c", 2),
     ("uo24", "Urine output in preceding 24 h, mL", "c", 0), ("imv", "Invasive mechanical ventilation", "b", 0),
     ("sedation_any", "Sedative infusion", "b", 0), ("hydrocort", "Intravenous hydrocortisone", "b", 0),
     ("dnr_by_stop", "Do-not-resuscitate order", "b", 0)]
rows = []
for lab, D, T, P in [("MIMIC-IV", d, t, pr), ("eICU-CRD", e, te, pre)]:
    mi = D.iloc[np.r_[P[:, 0], P[:, 1]]]; mt = np.r_[np.ones(len(P), int), np.zeros(len(P), int)]
    for c, name, kind, dec in V:
        if c not in D or D[c].notna().sum() == 0: continue
        b0, b1 = D[c][T == 0], D[c][T == 1]
        a0, a1 = pd.Series(mi[c].values[mt == 0]), pd.Series(mi[c].values[mt == 1])
        rows.append(dict(database=lab, variable=name, before_tapered=summ(b0, kind, dec), before_abrupt=summ(b1, kind, dec),
                         smd_before=round(smd(b1, b0, kind), 2), after_tapered=summ(a0, kind, dec),
                         after_abrupt=summ(a1, kind, dec), smd_after=round(smd(a1, a0, kind), 2)))
S = pd.DataFrame(rows); S.to_csv(RESULTS / "tableS4_psm_balance.csv", index=False)
print(S[~S.variable.str.contains("defines exposure")].groupby("database")[["smd_before", "smd_after"]].agg(lambda x: x.abs().max()))

# Missing data by discontinuation strategy (Table S2 footnote)
miss = []
for lab, D in (("MIMIC-IV", d), ("eICU-CRD", e)):
    for c in ["lac_last", "lac_max", "ph_last", "uo24", "apache", "frac_map65_24"]:
        if c in D and D[c].notna().any():
            miss.append(dict(database=lab, variable=c,
                             missing_tapered_pct=round(100 * D[c][D.abrupt == 0].isna().mean(), 1),
                             missing_abrupt_pct=round(100 * D[c][D.abrupt == 1].isna().mean(), 1)))
pd.DataFrame(miss).to_csv(RESULTS / "tableS2_missing_by_group.csv", index=False)
print(pd.DataFrame(miss))
