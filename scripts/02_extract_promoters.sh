#!/bin/bash -l
#SBATCH -p intel
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem 64gb
#SBATCH --time 02:00:00
#SBATCH --job-name=prom_extract
#SBATCH --out logs/prom_extract.%j.out
#SBATCH --err logs/prom_extract.%j.err

set -euo pipefail
source ~/.bashrc

SCRIPT="/bigdata/stajichlab/sadikshs/Aspergillus/scripts/extract_promoters.py"

# Your panel list WITHOUT the duplicate Af293
PANEL="/bigdata/stajichlab/sadikshs/Aspergillus/Asper_promoter_rh/promoter_rh/species_panel_step1/species_panel.no_ref.txt"

DNA_DIR="/bigdata/stajichlab/sadikshs/Aspergillus/genome_annotation/DNA"
GFF_DIR="/bigdata/stajichlab/sadikshs/Aspergillus/genome_annotation/gff"

OUT500="/bigdata/stajichlab/sadikshs/Aspergillus/promoters_new/by_species_up500"
OUT800="/bigdata/stajichlab/sadikshs/Aspergillus/promoters_new/by_species_up800"
mkdir -p logs "$OUT500" "$OUT800"

ID_ATTR="protein_id"
MINLEN=50

# ---- reference (hard-coded) ----
REF_BASE="FungiDB-68_AfumigatusAf293_AnnotatedProteins"
REF_DNA="${DNA_DIR}/${REF_BASE}.scaffolds.fa"
REF_GFF="${GFF_DIR}/${REF_BASE}.gff3"

echo "[ref] ${REF_BASE}"

python "$SCRIPT" --genome "$REF_DNA" --gff "$REF_GFF" \
  --out "${OUT800}/${REF_BASE}.up800.fa" --upstream 800 \
  --species "AfumigatusAf293_AnnotatedProteins" --id_attr "$ID_ATTR" \
  --drop_ambig --min_len "$MINLEN"

python "$SCRIPT" --genome "$REF_DNA" --gff "$REF_GFF" \
  --out "${OUT500}/${REF_BASE}.up500.fa" --upstream 500 \
  --species "AfumigatusAf293_AnnotatedProteins" --id_attr "$ID_ATTR" \
  --drop_ambig --min_len "$MINLEN"

# ---- helpers to find panel files by basename ----
find_dna () { ls -1 "${DNA_DIR}/${1}"*.fa* 2>/dev/null | head -n 1; }
find_gff () { ls -1 "${GFF_DIR}/${1}"*.gff* 2>/dev/null | head -n 1; }

while read -r proteome; do
  [[ -z "$proteome" ]] && continue
  [[ "${proteome:0:1}" == "#" ]] && continue

  base="${proteome%.proteins}"

  dna=$(find_dna "$base" || true)
  gff=$(find_gff "$base" || true)

  [[ -z "$dna" ]] && { echo "[WARN] missing DNA for $base (skip)" >&2; continue; }
  [[ -z "$gff" ]] && { echo "[WARN] missing GFF for $base (skip)" >&2; continue; }

  echo "[panel] $base"

  python "$SCRIPT" --genome "$dna" --gff "$gff" \
    --out "${OUT800}/${base}.up800.fa" --upstream 800 \
    --species "$base" --id_attr "$ID_ATTR" \
    --drop_ambig --min_len "$MINLEN"

  python "$SCRIPT" --genome "$dna" --gff "$gff" \
    --out "${OUT500}/${base}.up500.fa" --upstream 500 \
    --species "$base" --id_attr "$ID_ATTR" \
    --drop_ambig --min_len "$MINLEN"

done < "$PANEL"

echo "[done] by-species promoters:"
echo "  $OUT800"
echo "  $OUT500"