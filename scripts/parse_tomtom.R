#!/usr/bin/env Rscript
# parse_tomtom.R
# --------------
# Parses TomTom results from JASPAR and Cis-BP searches and produces:
#   1. matched_strict.tsv   — hits at E-value < 4e-3 
#   2. matched_relaxed.tsv  — hits at q-value <= 0.05 (for heatmap TF labels)
#   3. all_best_hits.tsv    — best hit per query regardless of significance
#   4. summary.txt          — match counts 
#
# Usage:
#   Rscript parse_tomtom.R <jaspar_tsv> <outdir> [cisbp_tsv]

suppressWarnings(suppressMessages({
  library(readr)
  library(dplyr)
  library(stringr)
}))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript parse_tomtom.R jaspar_tomtom.tsv outdir [cisbp_tomtom.tsv]",
       call. = FALSE)
}

jaspar_tsv <- args[1]
outdir     <- args[2]
cisbp_tsv  <- if (length(args) >= 3) args[3] else NULL

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ── Thresholds ────────────────────────────────
EVALUE_STRICT  <- 4e-3    # for counting matched motifs 
QVALUE_RELAXED <- 0.05    # for heatmap TF name labels 

# ── Helper to load and clean a TomTom TSV ──────────────────────
load_tomtom <- function(path, db_label) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  tt <- read_tsv(path, comment = "#", show_col_types = FALSE)
  names(tt) <- str_replace_all(names(tt), "[ -]", "_")
  tt <- tt %>%
    filter(!is.na(Query_ID), !is.na(Target_ID)) %>%
    mutate(
      database    = db_label,
      target_acc  = str_extract(Target_ID, "^[A-Z0-9]+\\.\\d+"),
      target_name = str_trim(str_remove(Target_ID, "^[A-Z0-9]+\\.\\d+\\s*"))
    )
  cat(sprintf("  Loaded %d rows from %s (%s)\n", nrow(tt), db_label, path))
  tt
}

cat("Loading TomTom results...\n")
jaspar <- load_tomtom(jaspar_tsv, "JASPAR")
cisbp  <- load_tomtom(cisbp_tsv,  "CisBP")

# ── Combine databases ───────────────────────────────────────────
# Union: for each query, keep the best hit across both databases.

if (!is.null(jaspar) && !is.null(cisbp)) {
  combined <- bind_rows(jaspar, cisbp)
  cat("  Combined rows:", nrow(combined), "\n")
} else if (!is.null(jaspar)) {
  combined <- jaspar
  cat("  [WARN] Only JASPAR results available. Add Cis-BP for complete analysis.\n")
} else {
  stop("No valid TomTom results found.")
}

# ── Best hit per query (no threshold — for annotation purposes) ─
all_best <- combined %>%
  arrange(Query_ID, E_value, q_value) %>%
  group_by(Query_ID) %>%
  slice(1) %>%
  ungroup()

# ── Strict threshold: E-value < 4e-3 ───────────

matched_strict <- combined %>%
  filter(E_value < EVALUE_STRICT) %>%
  arrange(Query_ID, E_value) %>%
  group_by(Query_ID) %>%
  slice(1) %>%      # best hit per query across both databases
  ungroup()

# ── Relaxed threshold: q-value <= 0.05 (for heatmap labels) ────
# Used to assign TF names in clustered heatmaps .
matched_relaxed <- combined %>%
  filter(q_value <= QVALUE_RELAXED) %>%
  arrange(Query_ID, E_value) %>%
  group_by(Query_ID) %>%
  slice(1) %>%
  ungroup()

# Per-database counts at strict threshold 
strict_per_db <- combined %>%
  filter(E_value < EVALUE_STRICT) %>%
  group_by(database) %>%
  summarise(
    unique_queries = n_distinct(Query_ID),
    total_hits     = n(),
    .groups = "drop")

# ── Write outputs ───────────────────────────────────────────────
write_tsv(all_best,         file.path(outdir, "all_best_hits.tsv"))
write_tsv(matched_strict,   file.path(outdir, "matched_strict_E4e3.tsv"))
write_tsv(matched_relaxed,  file.path(outdir, "matched_relaxed_q05.tsv"))

# Summary report
summary_lines <- c(
  paste0("Total queries (motifs): ", n_distinct(combined$Query_ID)),
  paste0(""),
  paste0("--- Strict threshold (E-value < 4e-3, paper Fig 2b) ---"),
  paste0("Matched motifs (union of databases): ", nrow(matched_strict)),
  paste0(""),
  paste0("  Per database:"),
  paste(capture.output(print(as.data.frame(strict_per_db), row.names=FALSE)),
        collapse="\n"),
  paste0(""),
  paste0("--- Relaxed threshold (q-value <= 0.05, for heatmap labels) ---"),
  paste0("Matched motifs: ", nrow(matched_relaxed)),
  paste0(""),
  paste0("Top strict matches:"),
  paste(capture.output(
    print(head(matched_strict %>%
                 select(Query_ID, Target_ID, database, E_value, q_value, Overlap), 15),
          row.names=FALSE)),
    collapse="\n")
)

writeLines(summary_lines, file.path(outdir, "summary.txt"))
cat(paste(summary_lines, collapse="\n"), "\n")
cat("\nOutputs written to:", outdir, "\n")