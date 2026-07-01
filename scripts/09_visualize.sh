#!/bin/bash -l
#SBATCH -p intel
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=32gb
#SBATCH --time=03:00:00
#SBATCH --job-name=visualize
#SBATCH --out logs/09_visualize.%j.out
#SBATCH --err logs/09_visualize.%j.err

set -euo pipefail
source ~/bigdata/Aspergillus/Asper_promoter_rh/Aspergillus_RHIEPA/scripts/00_config.sh

mkdir -p logs
module purge
module load R/4.5.0

# Install required R packages if missing
Rscript -e '
pkgs <- c("uwot","igraph","FNN","ggplot2","dplyr","readr",
          "stringr","tidyr","RColorBrewer","patchwork")
for (p in pkgs) {
  if (!requireNamespace(p, quietly=TRUE)) {
    install.packages(p, repos="https://cran.r-project.org", quiet=TRUE)
    cat("Installed:", p, "\n")
  }
}
cat("All packages ready\n")
'

# Input files
GAF="/bigdata/stajichlab/shared/projects/MosesLab_projects/Aspergillus/Af293/FungiDB-68_AfumigatusAf293_GO.gaf.gz"
SLIM_OBO="/srv/projects/db/GO/current/goslim_aspergillus.obo"
BASIC_OBO="/srv/projects/db/GO/current/go-basic.obo"

echo "=== Checking input files ==="
ls "${RESULTS_DIR}/motifs.meme"      && echo "OK: motifs.meme"      || echo "MISSING: motifs.meme"
ls "${TOMTOM_JASPAR_DIR}/tomtom.tsv" && echo "OK: tomtom.tsv"       || echo "MISSING: tomtom.tsv"
ls "${RESULTS_DIR}/rep.csv"          && echo "OK: rep.csv"           || echo "MISSING: rep.csv"
ls "$OG_MAP_FILTERED"                && echo "OK: og_map_filtered"   || echo "MISSING: og_map"
ls "$GAF"                            && echo "OK: GAF"               || echo "MISSING: GAF"
ls "$SLIM_OBO"                       && echo "OK: goslim_aspergillus.obo" || echo "MISSING: slim OBO"
ls "$BASIC_OBO"                      && echo "OK: go-basic.obo"      || echo "MISSING: basic OBO"

echo ""
echo "=== Running visualization ==="

Rscript "${SCRIPTS_DIR}/visualize.R" \
    "${RESULTS_DIR}/motifs.meme" \
    "${TOMTOM_JASPAR_DIR}/tomtom.tsv" \
    "${RESULTS_DIR}/rep.csv" \
    "${RESULTS_DIR}/visualization" \
    "$OG_MAP_FILTERED" \
    "$GAF" \
    "$SLIM_OBO" \
    "$BASIC_OBO"

echo ""
echo "=== Done. Results in: ${RESULTS_DIR}/visualization ==="
ls "${RESULTS_DIR}/visualization/"
echo ""
echo "Individual UMAPs per GO slim:"
ls "${RESULTS_DIR}/visualization/umap_per_go_slim/" | wc -l