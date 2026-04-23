#!/usr/bin/env bash
# ipsae_strc_reassess.sh
# Runs Dunbrack 2025 ipSAE (github.com/DunbrackLab/IPSAE) on STRC AF3 top-ranked
# models to get binder-confidence metric robust to ipTM's disordered-region
# false-positive failure mode. Parameters PAE cutoff 10 Å, dist cutoff 10 Å
# per IPSAE/README.md AF3 example. Outputs to research/strc/models/artifacts/ipsae/.
#
# Lit rationale: Dunbrack 2025 (bioRxiv 2025.02.10.637595) shows ipSAE is single
# best predictor of binding vs ipTM/iPAE/actifpTM/pDockQ, "correlated well with
# kD" (paper abstract). AF3 absolute Kd calibration NOT feasible per ACS 2024
# (jcim.4c00976) and bioRxiv 2025.04.07.647682.

set -euo pipefail

IPSAE="$HOME/DeepResearch/tools/IPSAE/ipsae.py"
AF3_RESULTS="$HOME/DeepResearch/strc/af3-results"
OUT="$HOME/Brain/research/strc/models/artifacts/ipsae"
PAE_CUT=10
DIST_CUT=10

mkdir -p "$OUT"

JOBS=(
  job-ultramini-x-tmem145-gold        # h09 gate 3 primary; pruned GOLD interface
  job-ultramini-x-tmem145-full        # h09 gate 3 secondary; full TMEM145
  job-ultramini-homodimer             # #26 engineered homodimer baseline
  job8-nfatc1-calcineurin             # KNOWN binder positive control (CnA-CnB physiological)
  job-d-mini-strc-tectorin-zp         # KNOWN non-binder negative control
  job-a-mini-strc-piezo2              # KNOWN non-binder negative control
)

for dir in "${JOBS[@]}"; do
  d="$AF3_RESULTS/$dir"
  if [ ! -d "$d" ]; then echo "MISS: $dir"; continue; fi
  json=$(ls "$d"/*full_data_0.json 2>/dev/null | head -1)
  cif=$(ls "$d"/*model_0.cif 2>/dev/null | head -1)
  if [ -z "$json" ] || [ -z "$cif" ]; then echo "NO FILES: $dir"; continue; fi
  echo "=== $dir ==="
  (cd "$d" && python3.11 "$IPSAE" "$json" "$cif" $PAE_CUT $DIST_CUT)
  base=$(basename "$cif" .cif)
  for suf in "_${PAE_CUT}_${DIST_CUT}.txt" "_${PAE_CUT}_${DIST_CUT}_byres.txt" "_${PAE_CUT}_${DIST_CUT}.pml"; do
    if [ -f "$d/${base}${suf}" ]; then mv "$d/${base}${suf}" "$OUT/${dir}${suf}"; fi
  done
done

echo "Done. Output in $OUT"
