#!/bin/bash -l
#SBATCH -p intel
#SBATCH -N 1 -n 8
#SBATCH --mem 128gb
#SBATCH --time 03:00:00
#SBATCH --job-name=group_by_og
#SBATCH --out logs/03_group_by_og.%j.out
#SBATCH --err logs/03_group_by_og.%j.err
set -euo pipefail
source ~/bigdata/Aspergillus/Asper_promoter_rh/Aspergillus_RHIEPA/scripts/00_config.sh

conda activate "$ENV_PROMOTER"

echo "=== Grouping promoters by OG (up800 — for training) ==="
python "${SCRIPTS_DIR}/03_group_by_og.py" \
  --map     "$OG_MAP_FILTERED" \
  --in_dir  "$BY_SPECIES_UP800" \
  --out_dir "$OG_UP800" \
  --min_seqs 1 \
  --ext ".fa"

echo ""
echo "=== Grouping promoters by OG (up500 — for representations) ==="
python "${SCRIPTS_DIR}/03_group_by_og.py" \
  --map     "$OG_MAP_FILTERED" \
  --in_dir  "$BY_SPECIES_UP500" \
  --out_dir "$OG_UP500" \
  --min_seqs 1 \
  --ext ".fa"

conda deactivate

echo ""
echo "OG file counts:"
echo "  up800: $(ls "$OG_UP800"/*.fa | wc -l)"
echo "  up500: $(ls "$OG_UP500"/*.fa | wc -l)"
