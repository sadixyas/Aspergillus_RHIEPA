#!/bin/bash -l
#SBATCH -p gpu --gres=gpu:1 --cpus-per-task=8 --mem 64gb --out logs/make_promoter_rh_weights_layer_info.%j.log --time=7-00:00:00

set -euo pipefail

cd /bigdata/stajichlab/sadikshs/Aspergillus/Asper_promoter_rh/promoter_rh

source ~/.bashrc
conda activate alan_tf_gpu

#python promoter_reverse_homology.py promoters_500bp
python promoter_reverse_homology.py /bigdata/stajichlab/sadikshs/Aspergillus/og_promoters/up800

conda deactivate