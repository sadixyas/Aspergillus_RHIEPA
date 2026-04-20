#!/usr/bin/env bash
set -euo pipefail
# shellcheck disable=SC1091
source "$(dirname "$0")/00_config.sh"

# Unique log per run (timestamp + nanoseconds)
LOG="$LOG_DIR/setup_$(date +%Y%m%d_%H%M%S_%N).log"
exec > >(tee -a "$LOG") 2>&1

stamp(){ printf "[%s] %s\n" "$(date '+%F %T')" "$*"; }

ORTHO_TSV="$ORTHO_DIR/Orthogroups.tsv"

stamp "[1/5] Make species_all.txt from Orthogroups.tsv header"
head -1 "$ORTHO_TSV" | awk -F'\t' '{for(i=2;i<=NF;i++) print $i}' > "$SPECIES_ALL"
wc -l "$SPECIES_ALL" | awk '{print "species_all.txt lines:",$1}'

stamp "[2/5] Choose species set and build cols.csv"
if [[ "$SPECIES_MODE" == "15" ]]; then
  if [[ ! -s "$SPECIES_15_FILE" ]]; then
    echo "ERROR: SPECIES_MODE=15 but file not found or empty: $SPECIES_15_FILE" >&2
    exit 1
  fi
  cp -f "$SPECIES_15_FILE" "$SPECIES_USE"
  echo "# Using 15-species list -> $SPECIES_USE"
else
  cp -f "$SPECIES_ALL" "$SPECIES_USE"
  echo "# Using ALL species -> $SPECIES_USE"
fi
paste -sd, "$SPECIES_USE" > "$COLS_CSV"
echo "cols.csv chars: $(awk '{print length}' "$COLS_CSV")"

stamp "[3/5] Build clean 3-column OG–species–gene map (single gene only)"
# Output: OG<TAB>Species<TAB>GeneID  (exactly 3 columns)
awk -F'\t' '
  FNR==NR {                      # pass 1: read chosen species list (one per line)
    s[++m]=$0; next
  }
  FNR==1 {                       # pass 2: THIS file (Orthogroups.tsv) header row
    for (i=1; i<=NF; i++) h[$i]=i;
    next
  }
  {
    og=$1
    for (i=1; i<=m; i++) {
      col = h[s[i]]
      if (!col) continue         # species not in header -> skip
      v = $col
      if (v=="") continue
      split(v, a, /,[[:space:]]*/)   # take the first gene if multiple
      if (a[1]!="") printf("%s\t%s\t%s\n", og, s[i], a[1])
    }
  }
' "$SPECIES_USE" "$ORTHO_TSV" > "$OG_MAP_3COL"

echo "3-col map lines: $(wc -l < "$OG_MAP_3COL")"

# Validate 3 columns right here (fail fast if broken)
awk -F'\t' 'NF!=3{print "ERROR: OG_MAP_3COL is not 3 columns at line",NR,"(NF="NF")"; exit 1} END{print "OG_MAP_3COL OK (3 columns)"}' "$OG_MAP_3COL"

stamp "[4/5] Select OGs with >= $THRESH single-copy presences (among chosen species)"
# Count directly from the wide matrix using chosen species (COLS_CSV)
awk -F'\t' -v COLS="$(cat "$COLS_CSV")" -v TH="$THRESH" '
  BEGIN{ n=split(COLS,c,",") }
  NR==1{ next }  # skip header
  {
    sc=0
    for(i=1;i<=n;i++){
      v=$c[i]
      if(v!="" && index(v,",")==0) sc++   # count only single-copy cells
    }
    if(sc>=TH) print $1
  }
' "$ORTHO_TSV" | sort -u > "$OG_SUBSET"
echo "OGs passing threshold: $(wc -l < "$OG_SUBSET")"

stamp "[5/5] Filter the 3-col map to those OGs only (final MAP)"
awk -F'\t' 'NR==FNR{keep[$1]=1; next} ($1 in keep)' \
  "$OG_SUBSET" "$OG_MAP_3COL" > "$OG_MAP_FILTERED"
echo "Filtered map lines: $(wc -l < "$OG_MAP_FILTERED")"

# FINAL VALIDATION — THIS IS THE CRITICAL FIX or else the script fails because of id mismatch 
awk -F'\t' 'NF!=3{print "ERROR: OG_MAP_FILTERED is not 3 columns at line",NR,"(NF="NF")"; exit 1} END{print "OG_MAP_FILTERED OK (3 columns)"}' "$OG_MAP_FILTERED"

# Small preview to sanity-check IDs
echo "Preview:"
head -3 "$OG_MAP_FILTERED" | sed 's/\t/ | /g'

stamp "Done. Final MAP file: $OG_MAP_FILTERED"

# Run time
if [[ ! ${START_TS-} ]]; then START_TS=$(date +%s); fi
END_TS=$(date +%s); DUR=$((END_TS - START_TS))
stamp "Setup completed in ${DUR}s"
stamp "Setup log saved to $LOG"
