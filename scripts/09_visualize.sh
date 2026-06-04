#!/bin/bash -l
#SBATCH -p intel
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=32gb
#SBATCH --time=02:00:00
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
          "stringr","tidyr","RColorBrewer")
for (p in pkgs) {
  if (!requireNamespace(p, quietly=TRUE)) {
    install.packages(p, repos="https://cran.r-project.org", quiet=TRUE)
    cat("Installed:", p, "\n")
  }
}
cat("All packages ready\n")
'

GAF="$HOME/bigdata/Aspergillus/Asper_promoter_rh/Aspergillus_RHIEPA"
GAF="/bigdata/stajichlab/shared/projects/MosesLab_projects/Aspergillus/Af293/FungiDB-68_AfumigatusAf293_GO.gaf.gz"
OBO="$HOME/bigdata/Aspergillus/go-basic.obo"

echo "Input files:"
ls "${RESULTS_DIR}/motifs.meme"      && echo "OK: motifs.meme"      || echo "MISSING: motifs.meme"
ls "${TOMTOM_JASPAR_DIR}/tomtom.tsv" && echo "OK: tomtom.tsv"       || echo "MISSING: tomtom.tsv"
ls "${RESULTS_DIR}/rep.csv"          && echo "OK: rep.csv"           || echo "MISSING: rep.csv"
ls "$OG_MAP_FILTERED"                && echo "OK: og_map_filtered"   || echo "MISSING: og_map"
ls "$GAF"                            && echo "OK: GAF"               || echo "MISSING: GAF"
ls "$OBO"                            && echo "OK: OBO"               || echo "MISSING: OBO"

Rscript "${SCRIPTS_DIR}/visualize.R" \
    "${RESULTS_DIR}/motifs.meme" \
    "${TOMTOM_JASPAR_DIR}/tomtom.tsv" \
    "${RESULTS_DIR}/rep.csv" \
    "${RESULTS_DIR}/visualization" \
    "$OG_MAP_FILTERED" \
    "$GAF" \
    "$OBO"

echo "Done. Results in: ${RESULTS_DIR}/visualization"
ls "${RESULTS_DIR}/visualization/"