#!/bin/bash -l
#SBATCH -p batch
#SBATCH -N 1 -n 8
#SBATCH --mem 32gb
#SBATCH --time 04:00:00
#SBATCH --job-name=tomtom
#SBATCH --out logs/07_tomtom.%j.out
#SBATCH --err logs/07_tomtom.%j.err
set -euo pipefail
source ~/bigdata/Aspergillus/Asper_promoter_rh/Aspergillus_RHIEPA/scripts/00_config.sh

conda activate "$ENV_MEME"

QUERY="${RESULTS_DIR}/motifs.meme"
[[ -s "$QUERY" ]] || { echo "[ERROR] Missing $QUERY"; exit 1; }

# ── Run TomTom against JASPAR fungi ───────────────────────────
[[ -s "$DB_JASPAR" ]] || { echo "[ERROR] Missing JASPAR DB: $DB_JASPAR"; exit 1; }
echo "=== TomTom vs JASPAR fungi ==="
rm -rf "$TOMTOM_JASPAR_DIR"
mkdir -p "$TOMTOM_JASPAR_DIR"
time tomtom \
  -oc      "$TOMTOM_JASPAR_DIR" \
  -evalue \
  -thresh  "$TOMTOM_THRESH" \
  -min-overlap 5 \
  "$QUERY" "$DB_JASPAR"
echo "JASPAR hits: $(tail -n +2 "${TOMTOM_JASPAR_DIR}/tomtom.tsv" | grep -v "^#" | wc -l)"

# ── Run TomTom against Cis-BP 2.0 (Aspergillus) ───────────────
# Download from: https://cisbp.ccbr.utoronto.ca
#   Species: Aspergillus fumigatus (or closest available)
#   Format: MEME
# Save to: $DB_CISBP (set in 00_config.sh)
if [[ -s "$DB_CISBP" ]]; then
  echo ""
  echo "=== TomTom vs Cis-BP 2.0 Aspergillus ==="
  rm -rf "$TOMTOM_CISBP_DIR"
  mkdir -p "$TOMTOM_CISBP_DIR"
  time tomtom \
    -oc      "$TOMTOM_CISBP_DIR" \
    -evalue \
    -thresh  "$TOMTOM_THRESH" \
    -min-overlap 5 \
    "$QUERY" "$DB_CISBP"
  echo "Cis-BP hits: $(tail -n +2 "${TOMTOM_CISBP_DIR}/tomtom.tsv" | grep -v "^#" | wc -l)"
else
  echo "[WARN] Cis-BP database not found at: $DB_CISBP"
  echo "       Download from https://cisbp.ccbr.utoronto.ca"
  echo "       Aspergillus fumigatus or closest available species."
  echo "       Skipping Cis-BP search (results will be incomplete)."
fi

conda deactivate
echo ""
echo "TomTom complete. Results in:"
echo "  JASPAR: $TOMTOM_JASPAR_DIR"
echo "  Cis-BP: $TOMTOM_CISBP_DIR"