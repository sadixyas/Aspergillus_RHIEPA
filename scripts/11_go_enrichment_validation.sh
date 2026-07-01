#!/bin/bash -l
#SBATCH -p intel
#SBATCH -N 1 -n 4
#SBATCH --mem=32gb
#SBATCH --time=02:00:00
#SBATCH --job-name=tf_go_enrich
#SBATCH --out logs/11_go_enrichment.%j.out
#SBATCH --err logs/11_go_enrichment.%j.err

set -euo pipefail
source ~/bigdata/Aspergillus/Asper_promoter_rh/Aspergillus_RHIEPA/scripts/00_config.sh

mkdir -p logs
module purge
module load R/4.5.0

# Install required packages if missing
Rscript -e '
pkgs <- c("readr","dplyr","stringr","tidyr","ggplot2","purrr","RColorBrewer")
for (p in pkgs) {
  if (!requireNamespace(p, quietly=TRUE)) {
    install.packages(p, repos="https://cran.r-project.org", quiet=TRUE)
    cat("Installed:", p, "\n")
  }
}
cat("Packages ready\n")
'

REP_CSV="${RESULTS_DIR}/rep.csv"
UMAP_GO="${RESULTS_DIR}/visualization/umap_coordinates_with_go.tsv"
TOMTOM_STRICT="${PARSE_DIR}/matched_strict_E4e3.tsv"
IC_TBL="${RESULTS_DIR}/visualization/motif_information_content.tsv"
OUTDIR="${RESULTS_DIR}/tf_go_enrichment_validation"

echo "=== Checking input files ==="
ls "$REP_CSV"       && echo "OK: rep.csv"
ls "$UMAP_GO"        && echo "OK: umap_coordinates_with_go.tsv"
ls "$TOMTOM_STRICT"  && echo "OK: matched_strict_E4e3.tsv"
ls "$IC_TBL"         && echo "OK: motif_information_content.tsv"

echo ""
echo "=== Running TF x GO enrichment validation ==="

Rscript "${SCRIPTS_DIR}/go_enrichment_validation.R" \
    "$REP_CSV" \
    "$UMAP_GO" \
    "$TOMTOM_STRICT" \
    "$IC_TBL" \
    "$OUTDIR" \
    4.0 \
    50 \
    0.05

echo ""
echo "=== Done. Results in: $OUTDIR ==="
ls "$OUTDIR/"
echo ""
echo "--- Summary report ---"
cat "$OUTDIR/SUMMARY_REPORT.txt"