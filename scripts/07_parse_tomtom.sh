#!/bin/bash -l
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 8gb
#SBATCH --time 00:30:00
#SBATCH --out logs/02_parse_tomtom.%j.out
#SBATCH --err logs/02_parse_tomtom.%j.err

set -euo pipefail

cd /bigdata/stajichlab/sadikshs/Aspergillus/Asper_promoter_rh/promoter_rh
mkdir -p logs

module purge
module load R/4.5.0

IN="tomtom_out/tomtom.tsv"
OUTDIR="results_01_tomtom"

Rscript parse_tomtom.R "$IN" "$OUTDIR"