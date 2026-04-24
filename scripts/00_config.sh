#!/usr/bin/env bash
# =============================================================
# 00_config.sh — single source of truth for the entire pipeline
# Source this file at the top of every other script:
#   source "$(dirname "$0")/00_config.sh"
# =============================================================

# ── Conda environments ────────────────────────────────────────
ENV_PROMOTER="promoter_tools"   # has biopython, gffutils, etc.
ENV_TF_GPU="alan_tf_gpu"        # tensorflow + GPU
ENV_MEME="meme_suite"           # MEME / TomTom

# ── Root directories ──────────────────────────────────────────
ROOT="$HOME/bigdata/Aspergillus"
SCRIPTS_DIR="$ROOT/Asper_promoter_rh/scripts"   # location of this file

# ── Input data ────────────────────────────────────────────────
ORTHO_DIR="$ROOT/Orthogroups"
GFF_DIR="$ROOT/genome_annotation/gff"
DNA_DIR="$ROOT/genome_annotation/DNA"

# Reference genome (Af293)
REF_BASE="FungiDB-68_AfumigatusAf293_AnnotatedProteins"
REF_DNA="${DNA_DIR}/${REF_BASE}.scaffolds.fa"
REF_GFF="${GFF_DIR}/${REF_BASE}.gff3"
REF_SPECIES="AfumigatusAf293"

# Panel file: one proteome basename per line, no header, no ref
PANEL="$SCRIPTS_DIR/species_panel.no_ref.txt"

# ── Intermediate outputs ──────────────────────────────────────
# Per-species promoter fastas (written by 02_extract_promoters.sh)
SPECIES_PROM_DIR="$ROOT/promoters_new"
BY_SPECIES_UP500="${SPECIES_PROM_DIR}/by_species_up500"
BY_SPECIES_UP800="${SPECIES_PROM_DIR}/by_species_up800"

# Per-OG promoter fastas (written by 03_group_by_og.py)
OG_PROM_DIR="$ROOT/og_promoters"
OG_UP800="${OG_PROM_DIR}/up800"
OG_UP500="${OG_PROM_DIR}/up500"

# Ortholog filtering outputs (written by 01_setup.sh)
SPECIES_ALL="${ORTHO_DIR}/species_all.txt"
SPECIES_USE="${ORTHO_DIR}/species_use.txt"
OG_MAP_3COL="${ORTHO_DIR}/og_species_gene_map.tsv"
OG_SUBSET="${ORTHO_DIR}/og_subset.txt"
OG_MAP_FILTERED="${ORTHO_DIR}/og_map_filtered.tsv"
ORTHO_TSV="${ORTHO_DIR}/Orthogroups.tsv"

# ── Ortholog filtering parameters ─────────────────────────────
# With 406 species, require ~65% coverage = 264 species minimum.
# Do NOT require single-copy — allows gene families with paralogs.
THRESH=264

# ── Model / training parameters ───────────────────────────────
FILTER_NUM=256      # number of PWMs to learn
FILTER_LEN=15       # PWM width
FAM_SIZE=4          # homologs averaged per contrastive step
BATCH_SIZE=256
NUM_EPOCHS=100
TRAIN_SEQ_L=800     # sequence length used during training
REP_SEQ_L=500       # sequence length used for representation 
MIN_ORTHOLOGS=5     # minimum orthologs to compute a representation

# Derived model filename (used by downstream scripts)
MODEL_STEM="promoter_rh_interpretableN${FILTER_NUM}L${FILTER_LEN}_f${FAM_SIZE}_b${BATCH_SIZE}_e${NUM_EPOCHS}"
MODEL_WEIGHTS="${MODEL_STEM}_final_weights.h5"
MATRIX_LAYER="${MODEL_STEM}_matrix_layer.txt"

# ── Motif database paths ──────────────────────────────────────
MOTIF_DB_DIR="$ROOT/Asper_promoter_rh/motif_db"
DB_JASPAR="${MOTIF_DB_DIR}/JASPAR2026_CORE_fungi_non-redundant_pfms_meme.txt"
DB_CISBP="${MOTIF_DB_DIR}/cisbp2_aspergillus.meme"   # download from cisbp.ccbr.utoronto.ca

# TomTom thresholds 
TOMTOM_THRESH=0.05          # passed to tomtom -thresh (E-value mode)
EVALUE_STRICT="0.004"       # 4e-3: used in parse step to count matches
QVALUE_RELAXED="0.05"       # used for heatmap TF name labels

# ── Output directories ────────────────────────────────────────
RH_DIR="$ROOT/Asper_promoter_rh/promoter_rh"   # working dir for training
RESULTS_DIR="${RH_DIR}/results"
TOMTOM_JASPAR_DIR="${RESULTS_DIR}/tomtom_jaspar"
TOMTOM_CISBP_DIR="${RESULTS_DIR}/tomtom_cisbp"
PARSE_DIR="${RESULTS_DIR}/parsed_tomtom"
LOG_DIR="${ROOT}/logs"

# ── GFF extraction parameters ─────────────────────────────────
ID_ATTR="protein_id"
MINLEN=50

# ── Create directories that must exist before any job runs ────
mkdir -p "$LOG_DIR" "$RESULTS_DIR" \
         "$BY_SPECIES_UP500" "$BY_SPECIES_UP800" \
         "$OG_UP800" "$OG_UP500" \
         "$TOMTOM_JASPAR_DIR" "$TOMTOM_CISBP_DIR" \
         "$PARSE_DIR"