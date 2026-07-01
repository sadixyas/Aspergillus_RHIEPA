#!/bin/bash -l
#SBATCH -p intel
#SBATCH -N 1 -n 4
#SBATCH --mem=32gb
#SBATCH --time=01:00:00
#SBATCH --job-name=clustering_figs
#SBATCH --out logs/10_cluster_figures.%j.out
#SBATCH --err logs/10_cluster_figures.%j.err

set -euo pipefail
source ~/bigdata/Aspergillus/Asper_promoter_rh/Aspergillus_RHIEPA/scripts/00_config.sh

mkdir -p logs
module purge
module load R/4.5.0

# Install patchwork if needed
Rscript -e '
pkgs <- c("ggplot2","dplyr","readr","stringr","tidyr","patchwork","RColorBrewer")
for (p in pkgs) {
  if (!requireNamespace(p, quietly=TRUE))
    install.packages(p, repos="https://cran.r-project.org", quiet=TRUE)
}
cat("Packages ready\n")
'

UMAP_TSV="${RESULTS_DIR}/visualization/umap_coordinates_with_go.tsv"
OUTDIR="${RESULTS_DIR}/cluster_figures"

[[ -s "$UMAP_TSV" ]] || { echo "[ERROR] Missing: $UMAP_TSV"; exit 1; }

echo "Input: $UMAP_TSV"
echo "Output: $OUTDIR"

Rscript "${SCRIPTS_DIR}/cluster_figures.R" "$UMAP_TSV" "$OUTDIR"

echo ""
echo "Done. Figures in: $OUTDIR"
ls "$OUTDIR/"