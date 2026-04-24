#!/bin/bash -l
#SBATCH -p intel --cpus-per-task=16 --mem 256gb
#SBATCH --time 6:30:00
#SBATCH --job-name=compute_rep
#SBATCH --out logs/05_compute_rep.%j.out
#SBATCH --err logs/05_compute_rep.%j.err
set -euo pipefail
source "$(dirname "$0")/00_config.sh"

cd "$RH_DIR"
conda activate "$ENV_TF_GPU"

echo "Computing representations from: $OG_UP500"
echo "Weights:  $MODEL_WEIGHTS"
echo "SEQ_L (representation): $REP_SEQ_L"
echo "Min orthologs: $MIN_ORTHOLOGS"

python "${SCRIPTS_DIR}/compute_promoter_rep.py" \
  "$MODEL_WEIGHTS" \
  "$OG_UP500" \
  --seq_l      "$REP_SEQ_L" \
  --min_seqs   "$MIN_ORTHOLOGS" \
  --filter_num "$FILTER_NUM" \
  --filter_len "$FILTER_LEN" \
  --out        "${RESULTS_DIR}/rep.csv"

conda deactivate
echo "Done. Output: ${RESULTS_DIR}/rep.csv"