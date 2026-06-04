#!/bin/bash -l
#SBATCH -p batch
#SBATCH -N 1 -n 1
#SBATCH --mem 16gb
#SBATCH --time 01:00:00
#SBATCH --job-name=meme_convert
#SBATCH --out logs/06_meme.%j.out
#SBATCH --err logs/06_meme.%j.err
set -euo pipefail
source ~/bigdata/Aspergillus/Asper_promoter_rh/Aspergillus_RHIEPA/scripts/00_config.sh

cd "$RH_DIR"
module purge
module load R/4.5.0

IN="${RH_DIR}/${MATRIX_LAYER}"
OUT="${RESULTS_DIR}/motifs.meme"

[[ -s "$IN"  ]] || { echo "[ERROR] Missing matrix layer: $IN";  exit 1; }

Rscript "${SCRIPTS_DIR}/matrixlayer_to_meme.R" "$IN" "$OUT" "$FILTER_LEN" "$FILTER_NUM"

echo "DONE: wrote $OUT"
grep -c "^MOTIF" "$OUT" | xargs echo "  Motifs in MEME file:"