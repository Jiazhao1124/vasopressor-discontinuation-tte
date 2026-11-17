#!/usr/bin/env bash
# Reproduce all results and figures. Place the five query outputs in data/raw/ first.
set -euo pipefail
cd "$(dirname "$0")/analysis"
for s in 01_prepare_mimic 02_prepare_eicu 03_main_analysis 04_competing_risk \
         05_dose_response 06_descriptives 07_table1 08_psm_balance; do
  echo "== $s"; python "$s.py"
done
cd ../figures
mkdir -p output
python figureS1_flow.py
LANG=C.UTF-8 Rscript figure1_dose_response.R
LANG=C.UTF-8 Rscript figureS2_subgroups.R
