"""Standardised 24-h re-initiation risk by NEE at discontinuation (Figure 1).

Logistic model with a restricted cubic B-spline of the stop dose plus covariates,
g-computation over the observed covariate distribution, 200 bootstrap resamples.
Output: results/dose_response_mimic.csv, results/dose_response_eicu.csv
"""
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import SplineTransformer
from common import CONF_MIMIC, CONF_EICU, RESULTS, load_mimic, load_eicu

rng = np.random.default_rng(11)
d = load_mimic()
# The reported analysis drew 200 bootstrap index vectors for a preliminary competing-risk
# check from this generator before the dose-response bootstrap. Drawing them here keeps
# the random stream, and therefore the reported confidence bands, identical.
for _ in range(200):
    rng.integers(0, len(d), len(d))

GRID = [0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 0.1, 0.15, 0.2]


def dose_response(D, cols, catcol, y, n_boot=200, knots=5):
    X = D[cols].astype(float)
    X = X.fillna(X.median())
    X = pd.concat([X, pd.get_dummies(D[catcol]).astype(float)], axis=1).values
    X = (X - X.mean(0)) / (X.std(0) + 1e-9)
    xd = D[["last_d"]].values

    def fit(idx):
        st = SplineTransformer(n_knots=knots, degree=3, knots="quantile").fit(xd[idx])
        lr = LogisticRegression(max_iter=4000).fit(np.hstack([st.transform(xd[idx]), X[idx]]), y[idx])
        return [lr.predict_proba(np.hstack([st.transform(np.full((len(idx), 1), gv)), X[idx]]))[:, 1].mean()
                for gv in GRID]

    point = fit(np.arange(len(D)))
    bs = np.array([fit(rng.integers(0, len(D), len(D))) for _ in range(n_boot)])
    return pd.DataFrame({"end_dose": GRID, "risk": point,
                         "lo": np.percentile(bs, 2.5, 0), "hi": np.percentile(bs, 97.5, 0)})


m = dose_response(d, [c for c in CONF_MIMIC if c != "cu"], "cu",
                  (d.restart_gap_h <= 24).astype(int).values)
m.to_csv(RESULTS / "dose_response_mimic.csv", index=False)
print("MIMIC-IV\n", m.round(3).to_string(index=False))

e = load_eicu()
de = dose_response(e, [c for c in CONF_EICU if c != "ut"], "unittype",
                   (e.restart_gap_h <= 24).astype(int).values, knots=4)
de.to_csv(RESULTS / "dose_response_eicu.csv", index=False)
print("eICU-CRD\n", de.round(3).to_string(index=False))
