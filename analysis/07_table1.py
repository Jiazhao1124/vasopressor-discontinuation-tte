"""Table 1: baseline characteristics of the two cohorts (median (IQR) or n (%)).

Output: results/table1.csv, results/table1_missing.csv
"""
import pandas as pd
from common import RESULTS, load_mimic, load_eicu

d, e = load_mimic(), load_eicu()
CU_M = {"Cardiac Vascular Intensive Care Unit (CVICU)": "card", "Coronary Care Unit (CCU)": "card",
        "Medical Intensive Care Unit (MICU)": "med", "Medical/Surgical Intensive Care Unit (MICU/SICU)": "mix",
        "Surgical Intensive Care Unit (SICU)": "surg", "Trauma SICU (TSICU)": "surg"}
CU_E = {"Med-Surg ICU": "mix", "CCU-CTICU": "card", "Cardiac ICU": "card", "CTICU": "card",
        "CSICU": "card", "MICU": "med", "SICU": "surg"}
d["cu5"] = d.first_careunit.map(CU_M).fillna("oth")
e["cu5"] = e.unittype.map(CU_E).fillna("oth")
for D in (d, e):
    for k in ["card", "med", "surg", "mix", "oth"]:
        D["cu_" + k] = (D.cu5 == k).astype(int)


def med(s, dec=1):
    s = s.dropna(); f = f"{{:.{dec}f}}"
    return f"{f.format(s.median())} ({f.format(s.quantile(.25))}–{f.format(s.quantile(.75))})"


def cnt(s):
    s = s.dropna(); return f"{int(s.sum()):,} ({100*s.mean():.1f})"


rows = []
def H(t): rows.append((t, "", ""))
def R(lab, mc, ec, fn, dec=1, scale=1):
    f = (lambda D, c: med(D[c] * scale, dec)) if fn == "m" else (lambda D, c: cnt(D[c]))
    rows.append((lab, "—" if mc is None else f(d, mc), "—" if ec is None else f(e, ec)))

H("Demographics"); R("Age, years", "age", "age", "m"); R("Male sex", "male", "male", "c")
H("Illness severity"); R("SOFA score, ICU day 1", "sofa_d1", None, "m"); R("APACHE IVa score", None, "apache", "m")
R("APACHE IVa predicted hospital mortality, %", None, "apache_pred", "m", 1, 100); R("Sepsis-3", "sepsis3", None, "c")
H("ICU type"); R("Cardiac or cardiovascular", "cu_card", "cu_card", "c"); R("Medical", "cu_med", "cu_med", "c")
R("Surgical or trauma", "cu_surg", "cu_surg", "c"); R("Mixed medical–surgical", "cu_mix", "cu_mix", "c")
R("Neurological or other", "cu_oth", "cu_oth", "c")
H("Vasopressor support"); R("Episode duration, h", "ep_dur_h", "ep_dur_h", "m")
R("Time from ICU admission to discontinuation, h", "icu_h_at_stop", "icu_h_at_stop", "m")
R("Peak NEE, µg/kg/min", "max_d", "max_d", "m", 3); R("NEE at discontinuation, µg/kg/min", "last_d", "last_d", "m", 3)
R("Abrupt stop (NEE >0.05 µg/kg/min)", "abrupt", "abrupt", "c"); R("Vasopressin at discontinuation", "vaso_late", "vaso_late", "c")
H("Physiology in the 2 h before discontinuation"); R("Mean arterial pressure, mmHg", "map_m2", "map_m2", "m")
R("Heart rate, beats/min", "hr_m2", "hr_m2", "m"); R("Respiratory rate, breaths/min", "rr_m2", "rr_m2", "m")
R("SpO2, %", "spo2_m2", "spo2_m2", "m"); R("MAP readings <65 mmHg in preceding 24 h, %", "frac_map65_24", "frac_map65_24", "m", 1, 100)
H("Laboratory values, last before discontinuation"); R("Lactate, mmol/L", "lac_last", "lac_last", "m")
R("Peak lactate, mmol/L", "lac_max", "lac_max", "m"); R("Creatinine, mg/dL", "cr_last", "cr_last", "m")
R("Arterial pH", "ph_last", "ph_last", "m", 2); R("Urine output in preceding 24 h, mL", "uo24", None, "m", 0)
H("Treatments at discontinuation"); R("Invasive mechanical ventilation", "imv", "imv", "c")
R("Sedative infusion", "sedation_any", None, "c"); R("Intravenous hydrocortisone", "hydrocort", None, "c")
R("Do-not-resuscitate order", "dnr_by_stop", "dnr_by_stop", "c")

hdr = (f"MIMIC-IV (n = {len(d):,} stays)", f"eICU-CRD (n = {len(e):,} patients)")
pd.DataFrame(rows, columns=["characteristic", *hdr]).to_csv(RESULTS / "table1.csv", index=False)
miss = [(c, round(100 * d[c].isna().mean(), 1) if c in d else None,
         round(100 * e[c].isna().mean(), 1) if c in e else None)
        for c in ["lac_last", "lac_max", "ph_last", "apache", "frac_map65_24", "uo24"]]
pd.DataFrame(miss, columns=["variable", "missing_mimic_pct", "missing_eicu_pct"]).to_csv(RESULTS / "table1_missing.csv", index=False)
for r in rows: print(" | ".join(r))
print(miss)
