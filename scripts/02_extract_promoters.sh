#!/bin/bash -l
#SBATCH -p intel
#SBATCH -N 1 -n 8
#SBATCH --mem 64gb
#SBATCH --time 02:00:00
#SBATCH --job-name=prom_extract
#SBATCH --out logs/02_extract_promoters.%j.out
#SBATCH --err logs/02_extract_promoters.%j.err
set -euo pipefail
source ~/bigdata/Aspergillus/Asper_promoter_rh/Aspergillus_RHIEPA/scripts/00_config.sh

SCRIPT="${SCRIPTS_DIR}/extract_promoters.py"

find_dna() { ls -1 "${DNA_DIR}/${1}"*.fa* 2>/dev/null | head -n 1; }
find_gff() { ls -1 "${GFF_DIR}/${1}"*.gff* 2>/dev/null | head -n 1; }

# ── Reference species (Af293) ──────────────────────────────────
echo "[ref] ${REF_SPECIES}"
python "$SCRIPT" \
  --genome "$REF_DNA" --gff "$REF_GFF" \
  --out "${BY_SPECIES_UP800}/${REF_BASE}.up800.fa" --upstream 800 \
  --species "$REF_SPECIES" --id_attr "protein_source_id" \
  --drop_ambig --min_len "$MINLEN"

python "$SCRIPT" \
  --genome "$REF_DNA" --gff "$REF_GFF" \
  --out "${BY_SPECIES_UP500}/${REF_BASE}.up500.fa" --upstream 500 \
  --species "$REF_SPECIES" --id_attr "protein_source_id" \
  --drop_ambig --min_len "$MINLEN"

# ── Panel species ──────────────────────────────────────────────
while read -r proteome; do
  [[ -z "$proteome" || "${proteome:0:1}" == "#" ]] && continue
  base="${proteome%.proteins}"
  dna=$(find_dna "$base" || true)
  gff=$(find_gff "$base" || true)
  [[ -z "$dna" ]] && { echo "[WARN] missing DNA for $base (skip)" >&2; continue; }
  [[ -z "$gff" ]] && { echo "[WARN] missing GFF for $base (skip)" >&2; continue; }

  echo "[panel] $base"
  python "$SCRIPT" \
    --genome "$dna" --gff "$gff" \
    --out "${BY_SPECIES_UP800}/${base}.up800.fa" --upstream 800 \
    --species "$base" --id_attr "Parent" \
    --drop_ambig --min_len "$MINLEN"

  python "$SCRIPT" \
    --genome "$dna" --gff "$gff" \
    --out "${BY_SPECIES_UP500}/${base}.up500.fa" --upstream 500 \
    --species "$base" --id_attr "Parent" \
    --drop_ambig --min_len "$MINLEN"

done < "$PANEL"

echo "[done] Per-species promoters written to:"
echo "  $BY_SPECIES_UP800"
echo "  $BY_SPECIES_UP500"