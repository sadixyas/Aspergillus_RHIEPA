#!/bin/bash -l
#SBATCH -p intel
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 16gb
#SBATCH --time 02:00:00
#SBATCH --out logs/03_viz_early_run.%j.out
#SBATCH --err logs/03_viz_early_run.%j.err

set -euo pipefail
cd /bigdata/stajichlab/sadikshs/Aspergillus/Asper_promoter_rh/promoter_rh
mkdir -p logs

module purge
module load R/4.5.0

Rscript visualize.R \
  motifs.meme \
  tomtom_out/tomtom.tsv \
  rep.csv \
  results_early_viz