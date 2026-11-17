-- MIMIC-IV-Note v2.2: comfort care / withdrawal and hospice keywords in the discharge summary
-- Output: data/raw/mimic_note.csv
WITH v AS (SELECT stay_id, starttime, endtime, norepinephrine_equivalent_dose AS d FROM mimiciv_derived.norepinephrine_equivalent_dose WHERE norepinephrine_equivalent_dose > 0),
o AS (SELECT v.*, MAX(endtime) OVER (PARTITION BY stay_id ORDER BY starttime ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING) AS prev_end FROM v),
g AS (SELECT o.*, SUM(CASE WHEN prev_end IS NULL OR starttime > prev_end + interval '60 minutes' THEN 1 ELSE 0 END) OVER (PARTITION BY stay_id ORDER BY starttime) AS ep FROM o),
e AS (SELECT stay_id, ep, MIN(starttime) AS st, MAX(endtime) AS et FROM g GROUP BY 1,2),
e2 AS (SELECT e.*, LEAD(st) OVER (PARTITION BY stay_id ORDER BY ep) AS next_st FROM e),
f AS (SELECT e2.stay_id, e2.st, e2.et, e2.next_st, i.subject_id, i.hadm_id, i.intime, i.outtime, a.deathtime
  FROM e2 JOIN mimiciv_icu.icustays i USING (stay_id) JOIN mimiciv_hosp.admissions a ON a.hadm_id = i.hadm_id
  WHERE e2.ep = 1 AND e2.et - e2.st >= interval '6 hours' AND i.outtime > e2.et + interval '1 hour'
    AND (a.deathtime IS NULL OR a.deathtime > e2.et + interval '1 hour'))
SELECT f.stay_id, COUNT(n.note_id) AS n_notes,
  MAX(CASE WHEN n.text ~* '(comfort measures|comfort[- ]focused|comfort care|\mCMO\M|withdraw[a-z]* (of )?(care|support|life)|transitioned to comfort|terminal(ly)? extubat)' THEN 1 ELSE 0 END) AS cmo_note,
  MAX(CASE WHEN n.text ~* 'hospice' THEN 1 ELSE 0 END) AS hospice_note
FROM f LEFT JOIN mimiciv_note.discharge n ON n.hadm_id = f.hadm_id
GROUP BY f.stay_id;
