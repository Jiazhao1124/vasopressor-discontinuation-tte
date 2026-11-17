-- MIMIC-IV v3.0: charted code status (itemid 223758) relative to discontinuation, 28-day and hospital death
-- Output: data/raw/mimic_cs.csv
WITH v AS (SELECT stay_id, starttime, endtime, norepinephrine_equivalent_dose AS d FROM mimiciv_derived.norepinephrine_equivalent_dose WHERE norepinephrine_equivalent_dose > 0),
o AS (SELECT v.*, MAX(endtime) OVER (PARTITION BY stay_id ORDER BY starttime ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING) AS prev_end FROM v),
g AS (SELECT o.*, SUM(CASE WHEN prev_end IS NULL OR starttime > prev_end + interval '60 minutes' THEN 1 ELSE 0 END) OVER (PARTITION BY stay_id ORDER BY starttime) AS ep FROM o),
e AS (SELECT stay_id, ep, MIN(starttime) AS st, MAX(endtime) AS et FROM g GROUP BY 1,2),
e2 AS (SELECT e.*, LEAD(st) OVER (PARTITION BY stay_id ORDER BY ep) AS next_st FROM e),
f AS (SELECT e2.stay_id, e2.st, e2.et, e2.next_st, i.subject_id, i.hadm_id, i.intime, i.outtime, a.deathtime
  FROM e2 JOIN mimiciv_icu.icustays i USING (stay_id) JOIN mimiciv_hosp.admissions a ON a.hadm_id = i.hadm_id
  WHERE e2.ep = 1 AND e2.et - e2.st >= interval '6 hours' AND i.outtime > e2.et + interval '1 hour'
    AND (a.deathtime IS NULL OR a.deathtime > e2.et + interval '1 hour')),
cs AS (SELECT f.stay_id, ce.charttime, ce.value FROM f JOIN mimiciv_icu.chartevents ce ON ce.stay_id = f.stay_id AND ce.itemid = 223758)
SELECT f.stay_id,
  (SELECT value FROM cs WHERE cs.stay_id = f.stay_id AND cs.charttime <= f.et ORDER BY cs.charttime DESC LIMIT 1) AS code_at_stop,
  MAX(CASE WHEN cs.charttime <= f.et + interval '2 hours' AND cs.value ILIKE '%comfort%' THEN 1 ELSE 0 END) AS cmo_by_stop2h,
  MAX(CASE WHEN cs.charttime <= f.et + interval '24 hours' AND cs.value ILIKE '%comfort%' THEN 1 ELSE 0 END) AS cmo_by_24h,
  MAX(CASE WHEN cs.charttime <= f.et AND (cs.value ILIKE '%DNR%' OR cs.value ILIKE '%DNAR%') THEN 1 ELSE 0 END) AS dnr_by_stop,
  MAX(a.discharge_location) AS disch_loc,
  MAX(CASE WHEN COALESCE(a.deathtime, p.dod::timestamp) <= f.et + interval '28 days' THEN 1 ELSE 0 END) AS death28,
  MAX(a.hospital_expire_flag) AS hosp_death
FROM f LEFT JOIN cs ON cs.stay_id = f.stay_id
JOIN mimiciv_hosp.admissions a ON a.hadm_id = f.hadm_id
JOIN mimiciv_hosp.patients p ON p.subject_id = f.subject_id
GROUP BY f.stay_id, f.et;
