#!/bin/bash -l
#SBATCH -p intel
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem 96gb
#SBATCH --time 06:00:00
#SBATCH --out logs/01_tomtom.%j.out
#SBATCH --err logs/01_tomtom.%j.err

set -euo pipefail
set -x

cd /bigdata/stajichlab/sadikshs/Aspergillus/Asper_promoter_rh/promoter_rh
mkdir -p logs

module purge

source ~/.bashrc
conda activate meme_suite

echo "tomtom: $(command -v tomtom)"
tomtom -version || true

QUERY="motifs.meme"
DB="motif_db/JASPAR2026_CORE_fungi_non-redundant_pfms_meme.txt"
OUTDIR="tomtom_out"

[[ -s "$QUERY" ]] || { echo "[ERROR] Missing $QUERY"; exit 1; }
[[ -s "$DB" ]] || { echo "[ERROR] Missing $DB"; exit 1; }

rm -rf "$OUTDIR"
mkdir -p "$OUTDIR"

echo "Starting tomtom..."
time tomtom -oc "$OUTDIR" -evalue -thresh 0.05 -min-overlap 5 "$QUERY" "$DB"
echo "Finished tomtom."

ls -lh "$OUTDIR" | head -n 30
head -n 5 "$OUTDIR/tomtom.tsv"

conda deactivate