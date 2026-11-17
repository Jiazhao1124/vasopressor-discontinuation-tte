-- MIMIC-IV v3.0: additional covariates (sedation, laboratory, 24-h haemodynamics, care unit)
-- Output: data/raw/mimic_vx.csv
WITH v AS (
  SELECT stay_id, starttime, endtime, norepinephrine_equivalent_dose AS d
  FROM mimiciv_derived.norepinephrine_equivalent_dose WHERE norepinephrine_equivalent_dose > 0),
o AS (SELECT v.*, MAX(endtime) OVER (PARTITION BY stay_id ORDER BY starttime ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING) AS prev_end FROM v),
g AS (SELECT o.*, SUM(CASE WHEN prev_end IS NULL OR starttime > prev_end + interval '60 minutes' THEN 1 ELSE 0 END) OVER (PARTITION BY stay_id ORDER BY starttime) AS ep FROM o),
e AS (SELECT stay_id, ep, MIN(starttime) AS st, MAX(endtime) AS et FROM g GROUP BY 1,2),
e2 AS (SELECT e.*, LEAD(st) OVER (PARTITION BY stay_id ORDER BY ep) AS next_st FROM e),
f AS (
  SELECT e2.stay_id, e2.st, e2.et, e2.next_st, i.subject_id, i.hadm_id, i.intime, i.outtime, a.deathtime
  FROM e2 JOIN mimiciv_icu.icustays i USING (stay_id)
  JOIN mimiciv_hosp.admissions a ON a.hadm_id = i.hadm_id
  WHERE e2.ep = 1 AND e2.et - e2.st >= interval '6 hours'
    AND i.outtime > e2.et + interval '1 hour'
    AND (a.deathtime IS NULL OR a.deathtime > e2.et + interval '1 hour')),

sed AS (
  SELECT f.stay_id,
    MAX(CASE WHEN ie.itemid = 222168 THEN 1 ELSE 0 END) AS propofol_on,
    MAX(CASE WHEN ie.itemid IN (221744,225942) THEN 1 ELSE 0 END) AS fentanyl_on,
    MAX(CASE WHEN ie.itemid = 225150 THEN 1 ELSE 0 END) AS dexmed_on,
    MAX(CASE WHEN ie.itemid = 221668 THEN 1 ELSE 0 END) AS midaz_on
  FROM f JOIN mimiciv_icu.inputevents ie ON ie.stay_id = f.stay_id AND ie.starttime <= f.et AND ie.endtime >= f.et - interval '30 minutes'
   AND ie.itemid IN (222168,221744,225942,225150,221668) GROUP BY 1),
fl AS (
  SELECT f.stay_id,
    SUM(CASE WHEN ie.starttime > f.et - interval '6 hours' THEN ie.amount END) AS fluid_in6,
    SUM(ie.amount) AS fluid_in24
  FROM f JOIN mimiciv_icu.inputevents ie ON ie.stay_id = f.stay_id AND ie.starttime > f.et - interval '24 hours' AND ie.starttime <= f.et
   AND ie.amountuom = 'ml' GROUP BY 1),
uo24 AS (
  SELECT f.stay_id, SUM(u.urineoutput) AS uo24
  FROM f JOIN mimiciv_derived.urine_output u ON u.stay_id = f.stay_id AND u.charttime > f.et - interval '24 hours' AND u.charttime <= f.et GROUP BY 1),
hc AS (
  SELECT f.stay_id, MAX(1) AS hydrocort
  FROM f JOIN mimiciv_hosp.prescriptions rx ON rx.hadm_id = f.hadm_id AND LOWER(rx.drug) LIKE '%hydrocortisone%' AND rx.route IN ('IV','IV BOLUS','IVPCA','IV DRIP')
   AND rx.starttime <= f.et AND (rx.stoptime IS NULL OR rx.stoptime >= f.et) GROUP BY 1),
cr AS (
  SELECT DISTINCT ON (f.stay_id) f.stay_id, c.creatinine AS cr_last, c.bun AS bun_last
  FROM f JOIN mimiciv_derived.chemistry c ON c.hadm_id = f.hadm_id AND c.creatinine IS NOT NULL AND c.charttime <= f.et AND c.charttime > f.et - interval '36 hours'
  ORDER BY f.stay_id, c.charttime DESC),
cbc AS (
  SELECT DISTINCT ON (f.stay_id) f.stay_id, b.wbc AS wbc_last, b.hemoglobin AS hb_last, b.platelet AS plt_last
  FROM f JOIN mimiciv_derived.complete_blood_count b ON b.hadm_id = f.hadm_id AND b.charttime <= f.et AND b.charttime > f.et - interval '36 hours'
  ORDER BY f.stay_id, b.charttime DESC),
bga AS (
  SELECT DISTINCT ON (f.stay_id) f.stay_id, b.ph AS ph_last, b.baseexcess AS be_last, b.pao2fio2ratio AS pf_last
  FROM f JOIN mimiciv_derived.bg b ON b.hadm_id = f.hadm_id AND b.ph IS NOT NULL AND b.charttime <= f.et AND b.charttime > f.et - interval '24 hours'
  ORDER BY f.stay_id, b.charttime DESC),
vs24 AS (
  SELECT f.stay_id,
    AVG(x.mbp) AS map_24, MIN(x.mbp) AS mapmin_24,
    AVG(CASE WHEN x.charttime > f.et - interval '6 hours' THEN CASE WHEN x.mbp < 65 THEN 1.0 ELSE 0.0 END END) AS frac_map65_6,
    AVG(CASE WHEN x.charttime > f.et - interval '24 hours' THEN CASE WHEN x.mbp < 65 THEN 1.0 ELSE 0.0 END END) AS frac_map65_24,
    REGR_SLOPE(x.mbp, EXTRACT(epoch FROM x.charttime - f.et)/3600) AS map_slope24,
    AVG(x.heart_rate) AS hr_24, MAX(x.heart_rate) AS hrmax_24
  FROM f JOIN mimiciv_derived.vitalsign x ON x.stay_id = f.stay_id AND x.charttime > f.et - interval '24 hours' AND x.charttime <= f.et GROUP BY 1),
dt AS (
  SELECT f.stay_id,
    MAX(CASE WHEN v.starttime <= f.et - interval '12 hours' AND v.endtime > f.et - interval '12 hours' THEN v.d END) AS d_m12,
    SUM(CASE WHEN v.d <= 0.03 THEN EXTRACT(epoch FROM v.endtime - v.starttime)/3600 ELSE 0 END) AS h_lowdose,
    REGR_SLOPE(v.d, EXTRACT(epoch FROM v.starttime - f.et)/3600) FILTER (WHERE v.starttime > f.et - interval '12 hours') AS dose_slope12
  FROM f JOIN v ON v.stay_id = f.stay_id AND v.starttime >= f.st AND v.starttime < f.et GROUP BY 1)
SELECT f.stay_id, i.first_careunit,
  sed.propofol_on, sed.fentanyl_on, sed.dexmed_on, sed.midaz_on, fl.fluid_in6, fl.fluid_in24, uo24.uo24, COALESCE(hc.hydrocort,0) AS hydrocort,
  cr.cr_last, cr.bun_last, cbc.wbc_last, cbc.hb_last, cbc.plt_last, bga.ph_last, bga.be_last, bga.pf_last,
  vs24.map_24, vs24.mapmin_24, vs24.frac_map65_6, vs24.frac_map65_24, vs24.map_slope24, vs24.hr_24, vs24.hrmax_24,
  dt.d_m12, dt.h_lowdose, dt.dose_slope12
FROM f JOIN mimiciv_icu.icustays i ON i.stay_id = f.stay_id
LEFT JOIN sed ON sed.stay_id = f.stay_id LEFT JOIN fl ON fl.stay_id = f.stay_id LEFT JOIN uo24 ON uo24.stay_id = f.stay_id
LEFT JOIN hc ON hc.stay_id = f.stay_id LEFT JOIN cr ON cr.stay_id = f.stay_id LEFT JOIN cbc ON cbc.stay_id = f.stay_id
LEFT JOIN bga ON bga.stay_id = f.stay_id LEFT JOIN vs24 ON vs24.stay_id = f.stay_id LEFT JOIN dt ON dt.stay_id = f.stay_id
