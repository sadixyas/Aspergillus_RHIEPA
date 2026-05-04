#!/bin/bash -l
#SBATCH -p batch
#SBATCH -N 1 -n 1
#SBATCH --mem 8gb
#SBATCH --time 00:30:00
#SBATCH --job-name=parse_tomtom
#SBATCH --out logs/08_parse_tomtom.%j.out
#SBATCH --err logs/08_parse_tomtom.%j.err
set -euo pipefail
source "$(dirname "$0")/00_config.sh"

module purge
module load R/4.5.0

JASPAR_TSV="${TOMTOM_JASPAR_DIR}/tomtom.tsv"
CISBP_TSV="${TOMTOM_CISBP_DIR}/tomtom.tsv"

[[ -s "$JASPAR_TSV" ]] || { echo "[ERROR] Missing $JASPAR_TSV"; exit 1; }

# Pass Cis-BP results only if the file exists
if [[ -s "$CISBP_TSV" ]]; then
  Rscript "${SCRIPTS_DIR}/parse_tomtom.R" \
    "$JASPAR_TSV" "$PARSE_DIR" "$CISBP_TSV"
else
  echo "[WARN] No Cis-BP results found. Running with JASPAR only."
  Rscript "${SCRIPTS_DIR}/parse_tomtom.R" \
    "$JASPAR_TSV" "$PARSE_DIR"
fi

echo "Parse complete. Results in: $PARSE_DIR"
cat "${PARSE_DIR}/summary.txt"