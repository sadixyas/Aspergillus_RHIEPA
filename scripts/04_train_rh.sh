#!/bin/bash -l
#SBATCH -p exfab
#SBATCH -A exfab
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=6-00:00:00
#SBATCH --job-name=train_rh
#SBATCH --out logs/04_train_rh.%j.out
#SBATCH --err logs/04_train_rh.%j.err
set -euo pipefail
source ~/bigdata/Aspergillus/Asper_promoter_rh/Aspergillus_RHIEPA/scripts/00_config.sh

cd "$RH_DIR"
conda activate "$ENV_TF_GPU"

echo "Training on: $OG_UP800"
echo "Model stem:  $MODEL_STEM"
echo "SEQ_L (training): $TRAIN_SEQ_L"

python "${SCRIPTS_DIR}/promoter_reverse_homology.py" "$OG_UP800"

conda deactivate