#!/bin/bash -l
#SBATCH -p intel
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 96gb
#SBATCH --time 02:30:00
#SBATCH --out logs/00_matrixlayer_to_meme.%j.out
#SBATCH --err logs/00_matrixlayer_to_meme.%j.err

set -euo pipefail

cd /bigdata/stajichlab/sadikshs/Aspergillus/Asper_promoter_rh/promoter_rh
mkdir -p logs

module purge
module load R/4.5.0

IN="promoter_rh_interpretableN256L15_f4_b256_e100_matrix_layer.txt"
OUT="motifs.meme"
RSCRIPT="matrixlayer_to_meme.R"

[[ -s "$IN" ]] || { echo "[ERROR] Missing input: $IN"; exit 1; }
[[ -s "$RSCRIPT" ]] || { echo "[ERROR] Missing R script: $RSCRIPT"; exit 1; }

if [[ -s "$OUT" ]]; then
  echo "[SKIP] $OUT already exists. Remove to regenerate: rm -f $OUT"
  exit 0
fi

Rscript "$RSCRIPT" "$IN" "$OUT" 15 256

echo "DONE ✅ wrote $OUT"
head -n 20 "$OUT"
grep -n "^MOTIF" "$OUT" | head