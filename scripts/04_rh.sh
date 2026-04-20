#!/bin/bash -l
#SBATCH -p intel --cpus-per-task=16 --mem 256gb --out logs/make_rep.csv.%j.log --time 6:30:00


cd /bigdata/stajichlab/sadikshs/Aspergillus/Asper_promoter_rh/promoter_rh

source ~/.bashrc
conda activate alan_tf_gpu

python compute_promoter_rep.py promoter_rh_interpretableN256L15_f4_b256_e100_final_weights.h5 /bigdata/stajichlab/sadikshs/Aspergillus/og_promoters/up500/
#python compute_promoter_rep.py promoter_rh_interpretableN256L15_f4_b256_final_weights.h5 /bigdata/stajichlab/sadikshs/Aspergillus/promoters_800bp/

conda deactivate