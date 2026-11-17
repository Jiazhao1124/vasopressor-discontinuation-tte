-- eICU-CRD v2.0: first-unit vasopressor episode (>=6 h, gaps <=120 min), features, outcomes, care limitation
-- Output: data/raw/eicu_vs.csv
WITH pt AS (
  SELECT p.patientunitstayid, p.hospitalid, p.unittype, p.unitdischargeoffset, p.unitdischargestatus, p.hospitaldischargestatus,
         p.hospitaldischargeoffset, p.gender,
         CASE WHEN p.age = '> 89' THEN 90 WHEN p.age ~ '^[0-9]+$' THEN p.age::int END AS age,
         COALESCE(NULLIF(p.admissionweight,0), 80) AS wt
  FROM eicuii.patient p
  WHERE p.unitvisitnumber = 1 AND (CASE WHEN p.age = '> 89' THEN 90 WHEN p.age ~ '^[0-9]+$' THEN p.age::int END) >= 18),
inf0 AS (
  SELECT i.patientunitstayid, i.infusionoffset AS t, i.drugname, i.drugrate::numeric AS r, pt.wt,
    CASE
      WHEN i.drugname ~* '(norepinephrine|levophed)' AND i.drugname ~* 'mcg/kg/min' THEN i.drugrate::numeric
      WHEN i.drugname ~* '(norepinephrine|levophed)' AND i.drugname ~* 'mcg/min'    THEN i.drugrate::numeric / pt.wt
      WHEN i.drugname ~* 'epinephrine' AND i.drugname !~* 'nor' AND i.drugname ~* 'mcg/kg/min' THEN i.drugrate::numeric
      WHEN i.drugname ~* 'epinephrine' AND i.drugname !~* 'nor' AND i.drugname ~* 'mcg/min'    THEN i.drugrate::numeric / pt.wt
      WHEN i.drugname ~* '(phenylephrine|neo-?synephrine)' AND i.drugname ~* 'mcg/kg/min' THEN i.drugrate::numeric / 10
      WHEN i.drugname ~* '(phenylephrine|neo-?synephrine)' AND i.drugname ~* 'mcg/min'    THEN i.drugrate::numeric / pt.wt / 10
      WHEN i.drugname ~* 'dopamine' AND i.drugname ~* 'mcg/kg/min' THEN i.drugrate::numeric / 100
      WHEN i.drugname ~* 'vasopressin' AND i.drugname ~* 'units/min' THEN i.drugrate::numeric * 2.5
      WHEN i.drugname ~* 'vasopressin' AND i.drugname ~* 'units/hr'  THEN i.drugrate::numeric * 2.5 / 60
      WHEN i.drugname ~* 'ml/hr' AND i.drugamount::text ~ '^[0-9]*\.?[0-9]+$' AND i.volumeoffluid::text ~ '^[0-9]*\.?[0-9]+$' AND i.volumeoffluid::numeric > 0 THEN
        CASE
          WHEN i.drugname ~* '(norepinephrine|levophed)' THEN i.drugrate::numeric * i.drugamount::numeric / i.volumeoffluid::numeric * 1000 / 60 / pt.wt
          WHEN i.drugname ~* 'epinephrine' AND i.drugname !~* 'nor' THEN i.drugrate::numeric * i.drugamount::numeric / i.volumeoffluid::numeric * 1000 / 60 / pt.wt
          WHEN i.drugname ~* '(phenylephrine|neo-?synephrine)' THEN i.drugrate::numeric * i.drugamount::numeric / i.volumeoffluid::numeric * 1000 / 60 / pt.wt / 10
          WHEN i.drugname ~* 'dopamine' THEN i.drugrate::numeric * i.drugamount::numeric / i.volumeoffluid::numeric * 1000 / 60 / pt.wt / 100
          WHEN i.drugname ~* 'vasopressin' THEN i.drugrate::numeric * i.drugamount::numeric / i.volumeoffluid::numeric / 60 * 2.5
        END
    END AS nee_raw,
    CASE WHEN i.drugname ~* 'vasopressin' THEN 1 ELSE 0 END AS is_avp
  FROM eicuii.infusiondrug i JOIN pt USING (patientunitstayid)
  WHERE i.drugname ~* '(norepinephrine|levophed|epinephrine|vasopressin|dopamine|phenylephrine|neo-?synephrine)'
    AND i.drugrate ~ '^[0-9]*\.?[0-9]+$' AND i.infusionoffset >= 0),
inf AS (SELECT inf0.*, CASE WHEN nee_raw > 0 AND nee_raw <= 5 THEN nee_raw END AS nee FROM inf0),
pos AS (SELECT * FROM inf WHERE r > 0),
o AS (SELECT patientunitstayid, t, LAG(t) OVER (PARTITION BY patientunitstayid ORDER BY t) AS prev_t FROM (SELECT DISTINCT patientunitstayid, t FROM pos) z),
g AS (SELECT patientunitstayid, t, SUM(CASE WHEN prev_t IS NULL OR t - prev_t > 120 THEN 1 ELSE 0 END) OVER (PARTITION BY patientunitstayid ORDER BY t) AS ep FROM o),
e AS (SELECT patientunitstayid, ep, MIN(t) AS st, MAX(t) AS et, COUNT(*) AS npts FROM g GROUP BY 1,2),
e2 AS (SELECT e.*, LEAD(st) OVER (PARTITION BY patientunitstayid ORDER BY ep) AS next_st FROM e),
f AS (
  SELECT e2.*, pt.hospitalid, pt.unittype, pt.age, pt.gender, pt.wt, pt.unitdischargeoffset, pt.unitdischargestatus,
         pt.hospitaldischargestatus, pt.hospitaldischargeoffset
  FROM e2 JOIN pt USING (patientunitstayid)
  WHERE e2.ep = 1 AND e2.et - e2.st >= 360 AND e2.npts >= 3
    AND pt.unitdischargeoffset > e2.et + 60
    AND NOT (pt.unitdischargestatus = 'Expired' AND pt.unitdischargeoffset <= e2.et + 60)),
dose AS (
  SELECT f.patientunitstayid,
    MAX(i.nee) AS max_d,
    MAX(i.is_avp) AS vaso_any,
    MAX(CASE WHEN i.is_avp = 1 AND i.t >= f.et - 120 THEN 1 ELSE 0 END) AS vaso_late,
    COUNT(DISTINCT CASE WHEN i.t >= f.et - 120 THEN regexp_replace(lower(i.drugname),'\s*\(.*','') END) AS n_agents_late
  FROM f JOIN pos i ON i.patientunitstayid = f.patientunitstayid AND i.t BETWEEN f.st AND f.et GROUP BY 1),
lastdrug AS (   -- last positive NEE per drug within 90 min of stop, summed
  SELECT patientunitstayid, SUM(nee) AS last_d, COUNT(*) FILTER (WHERE nee IS NULL) AS last_unk
  FROM (SELECT DISTINCT ON (f.patientunitstayid, regexp_replace(lower(i.drugname),'\s*\(.*',''))
          f.patientunitstayid, i.nee
        FROM f JOIN pos i ON i.patientunitstayid = f.patientunitstayid AND i.t BETWEEN f.et - 90 AND f.et
        ORDER BY f.patientunitstayid, regexp_replace(lower(i.drugname),'\s*\(.*',''), i.t DESC) z GROUP BY 1),
d6 AS (
  SELECT patientunitstayid, SUM(nee) AS d_m6
  FROM (SELECT DISTINCT ON (f.patientunitstayid, regexp_replace(lower(i.drugname),'\s*\(.*',''))
          f.patientunitstayid, i.nee
        FROM f JOIN pos i ON i.patientunitstayid = f.patientunitstayid AND i.t BETWEEN f.et - 420 AND f.et - 300
        ORDER BY f.patientunitstayid, regexp_replace(lower(i.drugname),'\s*\(.*',''), i.t DESC) z GROUP BY 1),
vs AS (
  SELECT f.patientunitstayid,
    AVG(CASE WHEN v.observationoffset > f.et - 120 THEN v.systemicmean END) AS map_m2_art,
    MIN(CASE WHEN v.observationoffset > f.et - 120 THEN v.systemicmean END) AS mapmin_m2_art,
    AVG(CASE WHEN v.observationoffset > f.et - 120 THEN v.heartrate END) AS hr_m2,
    AVG(v.heartrate) AS hr_24, AVG(v.systemicmean) AS map_24_art,
    AVG(CASE WHEN v.systemicmean < 65 THEN 1.0 WHEN v.systemicmean IS NOT NULL THEN 0.0 END) AS frac_map65_24,
    AVG(CASE WHEN v.observationoffset > f.et - 120 THEN v.respiration END) AS rr_m2,
    AVG(CASE WHEN v.observationoffset > f.et - 120 THEN v.sao2 END) AS spo2_m2
  FROM f JOIN eicuii.vitalperiodic v ON v.patientunitstayid = f.patientunitstayid AND v.observationoffset > f.et - 1440 AND v.observationoffset <= f.et GROUP BY 1),
nibp AS (
  SELECT f.patientunitstayid,
    AVG(CASE WHEN a.observationoffset > f.et - 180 THEN a.noninvasivemean END) AS map_m3_nibp,
    AVG(a.noninvasivemean) AS map_24_nibp
  FROM f JOIN eicuii.vitalaperiodic a ON a.patientunitstayid = f.patientunitstayid AND a.observationoffset > f.et - 1440 AND a.observationoffset <= f.et GROUP BY 1),
lab AS (
  SELECT f.patientunitstayid,
    (ARRAY_AGG(l.labresult ORDER BY l.labresultoffset DESC) FILTER (WHERE l.labname = 'lactate'))[1] AS lac_last,
    MAX(l.labresult) FILTER (WHERE l.labname = 'lactate') AS lac_max,
    (ARRAY_AGG(l.labresult ORDER BY l.labresultoffset DESC) FILTER (WHERE l.labname = 'creatinine'))[1] AS cr_last,
    (ARRAY_AGG(l.labresult ORDER BY l.labresultoffset DESC) FILTER (WHERE l.labname = 'pH'))[1] AS ph_last
  FROM f JOIN eicuii.lab l ON l.patientunitstayid = f.patientunitstayid AND l.labresultoffset BETWEEN f.st - 360 AND f.et
   AND l.labname IN ('lactate','creatinine','pH') GROUP BY 1),
vent AS (
  SELECT f.patientunitstayid, 1 AS imv
  FROM f JOIN eicuii.respiratorycare r ON r.patientunitstayid = f.patientunitstayid
   AND r.ventstartoffset <= f.et AND (r.ventendoffset >= f.et OR r.ventendoffset = 0) AND r.ventstartoffset > 0
  GROUP BY 1),
cpl AS (
  SELECT f.patientunitstayid,
    MAX(CASE WHEN c.cplitemoffset <= f.et + 120 AND c.cplitemvalue ILIKE '%comfort%' THEN 1 ELSE 0 END) AS cmo_by_stop2h,
    MAX(CASE WHEN c.cplitemoffset <= f.et + 1440 AND c.cplitemvalue ILIKE '%comfort%' THEN 1 ELSE 0 END) AS cmo_by_24h,
    MAX(CASE WHEN c.cplitemoffset <= f.et AND c.cplitemvalue IN ('Do not resuscitate','No CPR','No intubation','No vasopressors/inotropes') THEN 1 ELSE 0 END) AS dnr_by_stop,
    MAX(CASE WHEN c.cplitemoffset <= f.et + 1440 AND c.cplitemvalue = 'No vasopressors/inotropes' THEN 1 ELSE 0 END) AS novaso_order
  FROM f JOIN eicuii.careplangeneral c ON c.patientunitstayid = f.patientunitstayid AND c.cplgroup = 'Care Limitation' GROUP BY 1),
ap AS (
  SELECT patientunitstayid, MAX(apachescore) AS apache, MAX(predictedhospitalmortality::numeric) AS apache_pred
  FROM eicuii.apachepatientresult WHERE apacheversion = 'IVa' AND predictedhospitalmortality ~ '^[0-9.]+$' GROUP BY 1)
SELECT f.patientunitstayid, f.hospitalid, f.unittype, f.age, CASE WHEN f.gender = 'Male' THEN 1 ELSE 0 END AS male,
  ap.apache, ap.apache_pred,
  f.et / 60.0 AS icu_h_at_stop, (f.et - f.st) / 60.0 AS ep_dur_h, f.npts,
  dose.max_d, d6.d_m6, lastdrug.last_d, lastdrug.last_unk, dose.vaso_any, dose.vaso_late, dose.n_agents_late,
  vs.map_m2_art, vs.mapmin_m2_art, vs.map_24_art, vs.frac_map65_24, vs.hr_m2, vs.hr_24, vs.rr_m2, vs.spo2_m2,
  nibp.map_m3_nibp, nibp.map_24_nibp,
  lab.lac_last, lab.lac_max, lab.cr_last, lab.ph_last, COALESCE(vent.imv, 0) AS imv,
  COALESCE(cpl.cmo_by_stop2h, 0) AS cmo_by_stop2h, COALESCE(cpl.cmo_by_24h, 0) AS cmo_by_24h,
  COALESCE(cpl.dnr_by_stop, 0) AS dnr_by_stop, COALESCE(cpl.novaso_order, 0) AS novaso_order,
  (f.next_st - f.et) / 60.0 AS restart_gap_h,
  CASE WHEN f.unitdischargestatus = 'Expired' THEN (f.unitdischargeoffset - f.et) / 60.0 END AS icu_death_h,
  (f.unitdischargeoffset - f.et) / 60.0 AS icu_out_h,
  CASE WHEN f.hospitaldischargestatus = 'Expired' THEN 1 ELSE 0 END AS hosp_death,
  CASE WHEN f.hospitaldischargestatus = 'Expired' THEN (f.hospitaldischargeoffset - f.et) / 60.0 END AS hosp_death_h
FROM f
LEFT JOIN dose USING (patientunitstayid) LEFT JOIN lastdrug USING (patientunitstayid) LEFT JOIN d6 USING (patientunitstayid)
LEFT JOIN vs USING (patientunitstayid) LEFT JOIN nibp USING (patientunitstayid) LEFT JOIN lab USING (patientunitstayid)
LEFT JOIN vent USING (patientunitstayid) LEFT JOIN cpl USING (patientunitstayid) LEFT JOIN ap USING (patientunitstayid)
