-- MIMIC-IV v3.0: first vasopressor episode (>=6 h), features at discontinuation (T0)
-- Output: data/raw/mimic_vs.csv (one row per ICU stay)
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
dose AS (
  SELECT f.stay_id, MAX(v.d) AS max_d,
    MAX(CASE WHEN v.starttime <= f.et - interval '6 hours' AND v.endtime > f.et - interval '6 hours' THEN v.d END) AS d_m6,
    MAX(CASE WHEN v.starttime <= f.et - interval '2 hours' AND v.endtime > f.et - interval '2 hours' THEN v.d END) AS d_m2,
    SUM(v.d * EXTRACT(epoch FROM (LEAST(v.endtime, f.et) - v.starttime)) / 60) AS dose_auc
  FROM f JOIN v ON v.stay_id = f.stay_id AND v.starttime >= f.st AND v.starttime < f.et GROUP BY 1),
lastd AS (
  SELECT DISTINCT ON (f.stay_id) f.stay_id, v.d AS last_d, EXTRACT(epoch FROM v.endtime - v.starttime)/60 AS last_seg_min
  FROM f JOIN v ON v.stay_id = f.stay_id AND v.starttime >= f.st AND v.starttime < f.et
  ORDER BY f.stay_id, v.endtime DESC),
ag AS (
  SELECT f.stay_id,
    MAX(CASE WHEN va.vasopressin > 0 THEN 1 ELSE 0 END) AS vaso_any,
    MAX(CASE WHEN va.vasopressin > 0 AND va.endtime >= f.et - interval '2 hours' THEN 1 ELSE 0 END) AS vaso_late,
    MAX(CASE WHEN COALESCE(va.epinephrine,0) > 0 OR COALESCE(va.phenylephrine,0) > 0 OR COALESCE(va.dopamine,0) > 0 THEN 1 ELSE 0 END) AS other_any,
    MAX(CASE WHEN (COALESCE(va.epinephrine,0) > 0 OR COALESCE(va.phenylephrine,0) > 0 OR COALESCE(va.dopamine,0) > 0) AND va.endtime >= f.et - interval '2 hours' THEN 1 ELSE 0 END) AS other_late,
    MAX(CASE WHEN COALESCE(va.dobutamine,0) > 0 OR COALESCE(va.milrinone,0) > 0 THEN 1 ELSE 0 END) AS inotrope_any
  FROM f JOIN mimiciv_derived.vasoactive_agent va ON va.stay_id = f.stay_id AND va.starttime < f.et AND va.endtime > f.st GROUP BY 1),
vs AS (
  SELECT f.stay_id,
    AVG(CASE WHEN x.charttime > f.et - interval '2 hours' THEN x.mbp END) AS map_m2,
    MIN(CASE WHEN x.charttime > f.et - interval '2 hours' THEN x.mbp END) AS mapmin_m2,
    AVG(CASE WHEN x.charttime <= f.et - interval '2 hours' THEN x.mbp END) AS map_m6_2,
    STDDEV(x.mbp) AS map_sd6,
    AVG(CASE WHEN x.charttime > f.et - interval '2 hours' THEN x.sbp END) AS sbp_m2,
    AVG(CASE WHEN x.charttime > f.et - interval '2 hours' THEN x.dbp END) AS dbp_m2,
    AVG(CASE WHEN x.charttime > f.et - interval '2 hours' THEN x.heart_rate END) AS hr_m2,
    AVG(CASE WHEN x.charttime <= f.et - interval '2 hours' THEN x.heart_rate END) AS hr_m6_2,
    AVG(CASE WHEN x.charttime > f.et - interval '2 hours' THEN x.resp_rate END) AS rr_m2,
    AVG(CASE WHEN x.charttime > f.et - interval '2 hours' THEN x.spo2 END) AS spo2_m2,
    MAX(x.temperature) AS tmax_m6
  FROM f JOIN mimiciv_derived.vitalsign x ON x.stay_id = f.stay_id AND x.charttime > f.et - interval '6 hours' AND x.charttime <= f.et GROUP BY 1),
lac AS (
  SELECT DISTINCT ON (f.stay_id) f.stay_id, b.lactate AS lac_last, EXTRACT(epoch FROM f.et - b.charttime)/3600 AS lac_age_h
  FROM f JOIN mimiciv_derived.bg b ON b.hadm_id = f.hadm_id AND b.lactate IS NOT NULL AND b.charttime <= f.et AND b.charttime > f.et - interval '24 hours'
  ORDER BY f.stay_id, b.charttime DESC),
lacmax AS (
  SELECT f.stay_id, MAX(b.lactate) AS lac_max
  FROM f JOIN mimiciv_derived.bg b ON b.hadm_id = f.hadm_id AND b.lactate IS NOT NULL AND b.charttime <= f.et AND b.charttime >= f.st - interval '6 hours' GROUP BY 1),
uo AS (
  SELECT f.stay_id, SUM(u.urineoutput) AS uo6
  FROM f JOIN mimiciv_derived.urine_output u ON u.stay_id = f.stay_id AND u.charttime > f.et - interval '6 hours' AND u.charttime <= f.et GROUP BY 1),
vent AS (
  SELECT f.stay_id, MAX(CASE WHEN ve.ventilation_status = 'InvasiveVent' THEN 1 ELSE 0 END) AS imv
  FROM f JOIN mimiciv_derived.ventilation ve ON ve.stay_id = f.stay_id AND ve.starttime <= f.et AND ve.endtime >= f.et GROUP BY 1),
s3 AS (SELECT stay_id, MAX(CASE WHEN sepsis3 THEN 1 ELSE 0 END) AS sepsis3 FROM mimiciv_derived.sepsis3 GROUP BY 1)
SELECT f.stay_id, f.subject_id, p.anchor_year_group,
  d.admission_age AS age, CASE WHEN d.gender = 'M' THEN 1 ELSE 0 END AS male,
  fs.sofa AS sofa_d1, COALESCE(s3.sepsis3, 0) AS sepsis3,
  EXTRACT(epoch FROM f.et - f.intime)/3600 AS icu_h_at_stop,
  EXTRACT(epoch FROM f.et - f.st)/3600 AS ep_dur_h,
  EXTRACT(hour FROM f.et) AS stop_hour,
  dose.max_d, dose.d_m6, dose.d_m2, dose.dose_auc, lastd.last_d, lastd.last_seg_min,
  ag.vaso_any, ag.vaso_late, ag.other_any, ag.other_late, ag.inotrope_any,
  vs.map_m2, vs.mapmin_m2, vs.map_m6_2, vs.map_sd6, vs.sbp_m2, vs.dbp_m2, vs.hr_m2, vs.hr_m6_2, vs.rr_m2, vs.spo2_m2, vs.tmax_m6,
  lac.lac_last, lac.lac_age_h, lacmax.lac_max, uo.uo6, COALESCE(vent.imv, 0) AS imv,
  EXTRACT(epoch FROM f.next_st - f.et)/3600 AS restart_gap_h,
  EXTRACT(epoch FROM f.deathtime - f.et)/3600 AS death_h,
  EXTRACT(epoch FROM f.outtime - f.et)/3600 AS icu_out_h,
  CASE WHEN f.next_st <= f.et + interval '24 hours' OR f.deathtime <= f.et + interval '24 hours' THEN 1 ELSE 0 END AS fail24
FROM f
JOIN mimiciv_hosp.patients p ON p.subject_id = f.subject_id
LEFT JOIN mimiciv_derived.icustay_detail d ON d.stay_id = f.stay_id
LEFT JOIN mimiciv_derived.first_day_sofa fs ON fs.stay_id = f.stay_id
LEFT JOIN s3 ON s3.stay_id = f.stay_id
LEFT JOIN dose ON dose.stay_id = f.stay_id
LEFT JOIN lastd ON lastd.stay_id = f.stay_id
LEFT JOIN ag ON ag.stay_id = f.stay_id
LEFT JOIN vs ON vs.stay_id = f.stay_id
LEFT JOIN lac ON lac.stay_id = f.stay_id
LEFT JOIN lacmax ON lacmax.stay_id = f.stay_id
LEFT JOIN uo ON uo.stay_id = f.stay_id
LEFT JOIN vent ON vent.stay_id = f.stay_id
