#!/usr/bin/env bash

# Also need to specify the biopython created in conda env
PYTHON_BIN="$HOME/.conda/envs/promoter_tools/bin/python"

# Paths 
ROOT=~/bigdata/Aspergillus
ORTHO_DIR="$ROOT/Orthogroups"
GFF_DIR="$ROOT/genome_annotation/gff"
DNA_DIR="$ROOT/genome_annotation/DNA"
PROM_DIR="$ROOT/promoters_500bp" #Uncomment for the full species run
#PROM_DIR="$ROOT/promoters_800bp" 

# Where these helper scripts live
SCRIPTS_DIR=~/bigdata/Aspergillus/Asper_promoter_rh/Aspergillus/scripts

# My existing Python extractor
PROM_PY=~/bigdata/Aspergillus/Asper_promoter_rh/Aspergillus/scripts/extract_promoters.py

# Logs go here
LOG_DIR="$ROOT/logs"
mkdir -p "$LOG_DIR"

# Parameters defined
# Choose species set: "15 for initial trial" or "all"
SPECIES_MODE="all"

# If SPECIES_MODE=15, provide a file listing those 15 species (one per line, EXACTLY matching Orthogroups.tsv header)
SPECIES_15_FILE="$SCRIPTS_DIR/species_15.txt"

# Threshold for “single-copy” presence across the selected species
THRESH=12

# Output files written by 01_setup.sh that will be used by 02_extract_promoters.sh
SPECIES_ALL="$ORTHO_DIR/species_all.txt"
SPECIES_USE="$ORTHO_DIR/species_use.txt"                 # resolved list actually used (15 or all)
COLS_CSV="$ORTHO_DIR/cols.csv"
OG_MAP_3COL="$ORTHO_DIR/og_species_gene_map.correct.tsv" # 3-column OG–species–gene map
OG_SUBSET="$ORTHO_DIR/og_subset.txt"                     # OGs with >= THRESH single-copy presences
OG_MAP_FILTERED="$ORTHO_DIR/og_species_gene_map.filtered.tsv"  # final input to Python