#!/usr/bin/env bash
#SBATCH -p intel
#SBATCH -N 1 -n 4
#SBATCH --mem 32gb
#SBATCH --time 01:00:00
#SBATCH --job-name=og_setup
#SBATCH --out logs/01_setup.%j.out
#SBATCH --err logs/01_setup.%j.err
set -euo pipefail
source ~/bigdata/Aspergillus/Asper_promoter_rh/Aspergillus_RHIEPA/scripts/00_config.sh

LOG="${LOG_DIR}/setup_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

stamp() { printf "[%s] %s\n" "$(date '+%F %T')" "$*"; }

# ── [1/4] Extract species list from Orthogroups.tsv header ────
stamp "1/4 Building species_all.txt"
head -1 "$ORTHO_TSV" \
  | awk -F'\t' '{for(i=2;i<=NF;i++) print $i}' \
  > "$SPECIES_ALL"
N_SPECIES=$(wc -l < "$SPECIES_ALL")
echo "  Total species: $N_SPECIES"
echo "  Coverage threshold (THRESH): $THRESH  (~$(echo "scale=0; $THRESH*100/$N_SPECIES" | bc)%)"
cp -f "$SPECIES_ALL" "$SPECIES_USE"

# ── [2/4] Build 3-column OG–species–gene map ──────────────────
# Takes first gene when multiple paralogs exist (as in paper).
# Does NOT require single-copy — paralogs are allowed.
stamp "2/4 Building OG–species–gene map (all presence, not single-copy only)"
awk -F'\t' '
  FNR==NR { s[++m]=$0; next }          # pass 1: chosen species list
  FNR==1  { for(i=1;i<=NF;i++) h[$i]=i; next }  # pass 2: header row
  {
    og=$1
    for (i=1;i<=m;i++) {
      col=h[s[i]]
      if (!col) continue
      v=$col
      if (v=="") continue
      split(v, a, /,[[:space:]]*/)      # take only the first gene if paralogs exist
      if (a[1]!="") printf("%s\t%s\t%s\n", og, s[i], a[1])
    }
  }
' "$SPECIES_USE" "$ORTHO_TSV" > "$OG_MAP_3COL"

N_MAP=$(wc -l < "$OG_MAP_3COL")
echo "  3-col map lines: $N_MAP"
awk -F'\t' 'NF!=3 {print "ERROR: not 3 cols at line",NR; exit 1}
            END    {print "  OG_MAP_3COL OK (3 columns)"}' "$OG_MAP_3COL"

# ── [3/4] Select OGs with >= THRESH species present ───────────
# KEY FIX: count any non-empty cell (was: single-copy only).
# This retains gene families with lineage-specific duplications,
# which are important for secondary metabolism, pathogenicity, etc.
stamp "3/4 Selecting OGs with >= $THRESH species present (any copy)"
awk -F'\t' -v COLS="$(paste -sd, "$SPECIES_USE")" -v TH="$THRESH" '
  BEGIN { n=split(COLS,c,",") }
  NR==1 { next }
  {
    present=0
    for (i=1;i<=n;i++) {
      if ($c[i]!="") present++          # any copy counts, not just single-copy
    }
    if (present>=TH) print $1
  }
' "$ORTHO_TSV" | sort -u > "$OG_SUBSET"

N_OG=$(wc -l < "$OG_SUBSET")
echo "  OGs passing threshold: $N_OG"

# Diagnostic: how many were lost to the old single-copy filter?
OLD_COUNT=$(awk -F'\t' -v COLS="$(paste -sd, "$SPECIES_USE")" -v TH="$THRESH" '
  BEGIN { n=split(COLS,c,",") }
  NR==1 { next }
  {
    sc=0
    for(i=1;i<=n;i++) {
      v=$c[i]
      if(v!="" && index(v,",")==0) sc++
    }
    if(sc>=TH) print $1
  }
' "$ORTHO_TSV" | wc -l)
echo "  OGs that would pass old single-copy filter: $OLD_COUNT"
echo "  OGs recovered by removing single-copy restriction: $((N_OG - OLD_COUNT))"

# ── [4/4] Filter the 3-col map to selected OGs ────────────────
stamp "4/4 Filtering map to selected OGs"
awk -F'\t' 'NR==FNR{keep[$1]=1; next} ($1 in keep)' \
  "$OG_SUBSET" "$OG_MAP_3COL" > "$OG_MAP_FILTERED"

N_FILT=$(wc -l < "$OG_MAP_FILTERED")
echo "  Filtered map lines: $N_FILT"
awk -F'\t' 'NF!=3 {print "ERROR: not 3 cols at line",NR; exit 1}
            END    {print "  OG_MAP_FILTERED OK"}' "$OG_MAP_FILTERED"

echo ""
echo "Preview:"
head -3 "$OG_MAP_FILTERED" | column -t -s $'\t'

stamp "Done. Log: $LOG"