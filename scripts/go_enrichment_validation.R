#!/usr/bin/env Rscript
# 11_tf_go_enrichment_validation.R
# ----------------------------------
# Systematic validation: for every high-confidence TF motif (TomTom match,
# IC-filtered), test whether its predicted target genes are enriched in
# each Aspergillus GO slim biological process category, using Fisher's
# exact test with FDR correction across all TF x category comparisons.
#
# This produces a publication-style supplementary table and figure
# answering: "do motifs matching known TFs predict targets enriched in
# the GO categories that TF is independently known to regulate?"
#
# Designed to be run AFTER:
#   - visualize.R              (produces umap_coordinates_with_go.tsv)
#   - predicted_TF_target_edges_FILTERED.tsv (IC + score filtered edges)
#   - motif_information_content.tsv
#
# Usage:
#   Rscript 11_tf_go_enrichment_validation.R \
#       rep.csv \
#       umap_coordinates_with_go.tsv \
#       matched_strict_E4e3.tsv \
#       motif_information_content.tsv \
#       outdir \
#       [ic_cutoff] [top_n_targets] [fdr_cutoff]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop(
    "Usage: Rscript 11_tf_go_enrichment_validation.R ",
    "rep.csv umap_coordinates_with_go.tsv matched_strict_E4e3.tsv ",
    "motif_information_content.tsv outdir [ic_cutoff=4.0] ",
    "[top_n_targets=50] [fdr_cutoff=0.05]\n",
    call. = FALSE
  )
}

rep_csv     <- args[1]
umap_tsv    <- args[2]
tomtom_tsv  <- args[3]
ic_tsv      <- args[4]
outdir      <- args[5]
IC_CUTOFF   <- if (length(args) >= 6) as.numeric(args[6]) else 4.0
TOP_N       <- if (length(args) >= 7) as.integer(args[7]) else 50
FDR_CUTOFF  <- if (length(args) >= 8) as.numeric(args[8]) else 0.05

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

suppressMessages({
  library(readr); library(dplyr); library(stringr)
  library(tidyr); library(ggplot2); library(purrr)
})

for (pkg in c("RColorBrewer")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cran.r-project.org", quiet = TRUE)
}

cat(strrep("=", 70), "\n")
cat("TF motif x GO slim category enrichment validation\n")
cat(strrep("=", 70), "\n\n")
cat(sprintf("Parameters: IC cutoff = %.1f bits | top N targets per TF = %d | FDR cutoff = %.2f\n\n",
            IC_CUTOFF, TOP_N, FDR_CUTOFF))

# =============================================================
# 1. Load inputs
# =============================================================
cat("Loading inputs...\n")

rep <- read_csv(rep_csv, show_col_types = FALSE)
og_ids <- str_remove(rep$FILE, "\\.fa$")
X <- rep %>% select(where(is.numeric)) %>% as.matrix()
X[is.na(X)] <- 0
rownames(X) <- og_ids
colnames(X) <- paste0("matrix", 0:(ncol(X) - 1))
cat(sprintf("  rep.csv: %d OGs x %d motifs\n", nrow(X), ncol(X)))

umap_go <- read_tsv(umap_tsv, show_col_types = FALSE)
cat(sprintf("  GO annotations: %d OGs, %d categories (incl. unannotated)\n",
            nrow(umap_go), length(unique(umap_go$slim_name))))

tomtom <- read_tsv(tomtom_tsv, show_col_types = FALSE) %>%
  rename_with(~ gsub("-", "_", .x))
cat(sprintf("  TomTom strict matches: %d\n", nrow(tomtom)))

ic_tbl <- read_tsv(ic_tsv, show_col_types = FALSE)
cat(sprintf("  Motif IC table: %d motifs\n", nrow(ic_tbl)))

# =============================================================
# 2. Build high-confidence TF motif set
#    (TomTom strict match AND IC above cutoff)
# =============================================================
cat("\nFiltering to high-confidence TF motifs...\n")

tf_motifs <- tomtom %>%
  left_join(ic_tbl, by = c("Query_ID" = "motif")) %>%
  filter(ic_total > IC_CUTOFF) %>%
  mutate(
    tf_name = str_extract(Target_ID, "^[^ ]+"),
    tf_name = str_remove(tf_name, "_3\\.00$")
  ) %>%
  distinct(Query_ID, tf_name, ic_total, E_value)

cat(sprintf("  High-confidence TF motifs (IC > %.1f): %d\n",
            IC_CUTOFF, nrow(tf_motifs)))

if (nrow(tf_motifs) == 0) {
  stop("No motifs passed the IC cutoff. Lower ic_cutoff and retry.")
}

write_tsv(tf_motifs, file.path(outdir, "high_confidence_tf_motifs.tsv"))

# =============================================================
# 3. Build predicted target gene sets (top N per motif)
# =============================================================
cat("\nBuilding predicted target gene sets (top", TOP_N, "OGs per motif)...\n")

get_top_targets <- function(motif_col) {
  if (!motif_col %in% colnames(X)) return(character(0))
  scores <- X[, motif_col]
  if (all(abs(scores) < 1e-6)) return(character(0))  # dead motif guard
  ord <- order(scores, decreasing = TRUE)
  rownames(X)[ord[seq_len(min(TOP_N, length(ord)))]]
}

target_sets <- map(tf_motifs$Query_ID, get_top_targets)
names(target_sets) <- tf_motifs$Query_ID

# Drop any motifs that turned out to be dead (all-zero scores)
keep <- sapply(target_sets, length) > 0
cat(sprintf("  Motifs with non-trivial target scores: %d / %d\n",
            sum(keep), length(keep)))
target_sets <- target_sets[keep]
tf_motifs   <- tf_motifs %>% filter(Query_ID %in% names(target_sets))

# =============================================================
# 4. Build GO category gene sets (background = all annotated OGs)
# =============================================================
cat("\nBuilding GO slim category gene sets...\n")

go_sets <- umap_go %>%
  filter(slim_name != "unannotated") %>%
  group_by(slim_name) %>%
  summarise(OGs = list(OG), n = n(), .groups = "drop") %>%
  filter(n >= 5)   # skip categories too small for a meaningful test

cat(sprintf("  GO categories with >=5 annotated OGs: %d\n", nrow(go_sets)))

universe <- og_ids   # all OGs are the statistical background

# =============================================================
# 5. Fisher's exact test for every TF motif x GO category pair
# =============================================================
cat("\nRunning Fisher's exact tests (TF motif x GO category)...\n")
cat(sprintf("  Total tests: %d motifs x %d categories = %d\n",
            length(target_sets), nrow(go_sets),
            length(target_sets) * nrow(go_sets)))

run_fisher <- function(target_og_set, go_og_set, universe) {
  in_target_in_go     <- length(intersect(target_og_set, go_og_set))
  in_target_not_go     <- length(target_og_set) - in_target_in_go
  not_target_in_go     <- length(go_og_set) - in_target_in_go
  not_target_not_go    <- length(universe) - length(union(target_og_set, go_og_set))

  tbl <- matrix(c(in_target_in_go, in_target_not_go,
                  not_target_in_go, not_target_not_go),
                nrow = 2)
  ft <- fisher.test(tbl, alternative = "greater")
  list(
    overlap         = in_target_in_go,
    target_size     = length(target_og_set),
    go_size         = length(go_og_set),
    expected_overlap = round(length(target_og_set) * length(go_og_set) / length(universe), 2),
    fold_enrichment = round(
      (in_target_in_go / length(target_og_set)) /
      (length(go_og_set) / length(universe)), 3),
    p_value         = ft$p.value
  )
}

results <- list()
idx <- 1
for (i in seq_along(target_sets)) {
  motif_id <- names(target_sets)[i]
  tf_row   <- tf_motifs %>% filter(Query_ID == motif_id) %>% slice(1)

  for (j in seq_len(nrow(go_sets))) {
    go_cat <- go_sets$slim_name[j]
    go_og  <- go_sets$OGs[[j]]

    res <- run_fisher(target_sets[[motif_id]], go_og, universe)

    results[[idx]] <- tibble(
      motif            = motif_id,
      tf_name          = tf_row$tf_name,
      tomtom_E_value   = tf_row$E_value,
      motif_IC         = round(tf_row$ic_total, 2),
      go_category      = go_cat,
      target_set_size  = res$target_size,
      go_category_size = res$go_size,
      overlap          = res$overlap,
      expected_overlap = res$expected_overlap,
      fold_enrichment  = res$fold_enrichment,
      p_value          = res$p_value
    )
    idx <- idx + 1
  }
}

enrichment_tbl <- bind_rows(results) %>%
  mutate(fdr = p.adjust(p_value, method = "BH")) %>%
  arrange(fdr, desc(fold_enrichment))

write_tsv(enrichment_tbl, file.path(outdir, "tf_go_enrichment_full_results.tsv"))
cat(sprintf("\n  Total TF x category tests completed: %d\n", nrow(enrichment_tbl)))
cat(sprintf("  Significant at FDR < %.2f: %d\n",
            FDR_CUTOFF, sum(enrichment_tbl$fdr < FDR_CUTOFF)))

# =============================================================
# 6. Significant results summary table (publication supplementary table)
# =============================================================
cat("\nBuilding significant results summary...\n")

sig_results <- enrichment_tbl %>%
  filter(fdr < FDR_CUTOFF, overlap >= 3, fold_enrichment > 1) %>%
  arrange(fdr)

write_tsv(sig_results, file.path(outdir, "tf_go_enrichment_SIGNIFICANT.tsv"))
cat(sprintf("  Significant, biologically meaningful hits (overlap>=3, FE>1): %d\n",
            nrow(sig_results)))

if (nrow(sig_results) > 0) {
  cat("\nTop 15 significant TF -> GO category enrichments:\n")
  print(sig_results %>%
          select(tf_name, go_category, overlap, target_set_size,
                go_category_size, fold_enrichment, fdr) %>%
          head(15),
        n = 15)
}

# =============================================================
# 7. Highlight specific literature-supported TFs (if present)
#    RIM101/PacC -> pathogenesis / cell wall
#    SFP1        -> ribosome biogenesis
# =============================================================
cat("\n--- Literature-anchored validation checks ---\n")

check_tf_category <- function(tf_pattern, go_pattern, label) {
  hit <- enrichment_tbl %>%
    filter(str_detect(tf_name, regex(tf_pattern, ignore_case = TRUE)),
           str_detect(go_category, regex(go_pattern, ignore_case = TRUE)))
  if (nrow(hit) == 0) {
    cat(sprintf("  [%s] No matching motif/category found in data.\n", label))
    return(NULL)
  }
  hit <- hit %>% arrange(fdr) %>% slice(1)
  cat(sprintf(
    "  [%s] TF=%s  GO=%s  overlap=%d/%d  fold-enrichment=%.2fx  FDR=%.3g\n",
    label, hit$tf_name, hit$go_category, hit$overlap, hit$target_set_size,
    hit$fold_enrichment, hit$fdr))
  hit
}

# Note: tf_name column contains JASPAR IDs (e.g. MA0368.1) not common names.
# MA0368.1 = RIM101 (PacC ortholog, pH response / pathogenesis)
# MA0378.2 = SFP1   (ribosome biogenesis / TORC1 signaling)
# MA0299.1 = GAL4   (carbon source / carbohydrate metabolism)
# MA2703.1 = PPR1   (pyrimidine metabolism)
# MA0289.1 = DAL80  (nitrogen catabolite repression)
rim101_check  <- check_tf_category("MA0368\\.1", "pathogenesis",        "RIM101 (MA0368.1) -> pathogenesis")
rim101_check2 <- check_tf_category("MA0368\\.1", "cell wall|morphology","RIM101 (MA0368.1) -> cell wall/morphology")
sfp1_check    <- check_tf_category("MA0378\\.2", "ribosome",            "SFP1 (MA0378.2) -> ribosome biogenesis")
gal4_check    <- check_tf_category("MA0299\\.1", "carbohydrate",        "GAL4 (MA0299.1) -> carbohydrate metabolism")
ppr1_check    <- check_tf_category("MA2703\\.1", "cellular amino acid|vitamin|metabol", "PPR1 (MA2703.1) -> amino acid/vitamin metabolism")
dal80_check   <- check_tf_category("MA0289\\.1", "transport|nitrogen|cellular amino",   "DAL80 (MA0289.1) -> nitrogen/transport")

anchored <- bind_rows(rim101_check, rim101_check2, sfp1_check,
                      gal4_check, ppr1_check, dal80_check)
if (nrow(anchored) > 0) {
  write_tsv(anchored, file.path(outdir, "literature_anchored_checks.tsv"))
}

# =============================================================
# 8. Publication figure: volcano-style plot of all tests
# =============================================================
cat("\nGenerating publication figure...\n")

enrichment_tbl <- enrichment_tbl %>%
  mutate(
    sig_label = case_when(
      fdr < FDR_CUTOFF & fold_enrichment > 1 ~ "Significant enrichment",
      TRUE ~ "Not significant"
    ),
    neg_log10_fdr = -log10(pmax(fdr, 1e-300))
  )

p_volcano <- ggplot(enrichment_tbl,
                    aes(x = fold_enrichment, y = neg_log10_fdr,
                        color = sig_label)) +
  geom_point(alpha = 0.5, size = 1.2) +
  geom_hline(yintercept = -log10(FDR_CUTOFF),
             linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("Significant enrichment" = "#d6604d",
                                "Not significant" = "grey70")) +
  scale_x_log10() +
  labs(
    title    = "TF motif -> GO slim category enrichment (all pairwise tests)",
    subtitle = sprintf(
      "%d motifs x %d GO categories = %d tests | %d significant at FDR < %.2f (BH-corrected)",
      length(unique(enrichment_tbl$motif)), nrow(go_sets),
      nrow(enrichment_tbl), sum(enrichment_tbl$fdr < FDR_CUTOFF), FDR_CUTOFF),
    x = "Fold enrichment (log scale)",
    y = expression(-log[10]~"(FDR-adjusted p-value)"),
    color = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    legend.position = "bottom"
  )

ggsave(file.path(outdir, "enrichment_volcano_plot.pdf"),
       p_volcano, width = 8, height = 6)
ggsave(file.path(outdir, "enrichment_volcano_plot.png"),
       p_volcano, width = 8, height = 6, dpi = 300)
cat("  Saved: enrichment_volcano_plot.pdf/.png\n")

# =============================================================
# 9. Publication figure: top significant hits bar chart
# =============================================================
if (nrow(sig_results) > 0) {
  top_n_plot <- sig_results %>% head(20) %>%
    mutate(label = sprintf("%s -> %s", tf_name, go_category),
           label = str_wrap(label, width = 40),
           label = factor(label, levels = rev(unique(label))))

  p_bar <- ggplot(top_n_plot, aes(x = fold_enrichment, y = label)) +
    geom_col(fill = "#2166ac", width = 0.7) +
    geom_text(aes(label = sprintf("overlap=%d/%d, FDR=%.2g",
                                  overlap, target_set_size, fdr)),
              hjust = -0.05, size = 2.6) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.6))) +
    labs(
      title    = "Top significant TF motif -> GO slim category enrichments",
      subtitle = sprintf("FDR < %.2f (Benjamini-Hochberg corrected, Fisher's exact test)",
                         FDR_CUTOFF),
      x        = "Fold enrichment",
      y        = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title    = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, color = "grey40"),
      panel.grid.major.y = element_blank()
    )

  ggsave(file.path(outdir, "top_significant_enrichments_barchart.pdf"),
         p_bar, width = 9, height = max(5, nrow(top_n_plot) * 0.35))
  ggsave(file.path(outdir, "top_significant_enrichments_barchart.png"),
         p_bar, width = 9, height = max(5, nrow(top_n_plot) * 0.35), dpi = 300)
  cat("  Saved: top_significant_enrichments_barchart.pdf/.png\n")
}

# =============================================================
# 10. Written summary report (for direct inclusion in lab notes / methods)
# =============================================================
cat("\nWriting summary report...\n")

report_lines <- c(
  "TF MOTIF x GO SLIM CATEGORY ENRICHMENT VALIDATION",
  strrep("=", 60),
  "",
  "METHODS",
  strrep("-", 60),
  "High-confidence TF motifs were defined as those with a TomTom",
  "match (E < 4e-3) to a known JASPAR/CIS-BP motif AND information",
  sprintf("content (IC) > %.1f bits, to exclude low-specificity motifs that", IC_CUTOFF),
  sprintf("match reference motifs spuriously (n = %d motifs retained).", length(unique(enrichment_tbl$motif))),
  "",
  sprintf("For each high-confidence motif, the top %d orthogroups (OGs) by", TOP_N),
  "RHIEPA motif score were defined as predicted regulatory targets.",
  "",
  "For each Aspergillus GO slim biological process category with >=5",
  "annotated OGs (from FungiDB GAF + go-basic.obo hierarchy traversal,",
  "category assigned by literature-informed priority ranking), a one-sided",
  "Fisher's exact test (alternative = 'greater') was performed to test",
  "whether predicted targets of each TF motif were enriched for that",
  "category, against a background of all orthogroups represented in",
  "rep.csv. P-values were corrected for multiple testing across all",
  "TF x category pairs using the Benjamini-Hochberg (FDR) procedure.",
  "",
  "RESULTS",
  strrep("-", 60),
  sprintf("Total TF motifs tested:          %d", length(unique(enrichment_tbl$motif))),
  sprintf("Total GO categories tested:      %d", nrow(go_sets)),
  sprintf("Total pairwise tests:            %d", nrow(enrichment_tbl)),
  sprintf("Significant at FDR < %.2f:       %d", FDR_CUTOFF, sum(enrichment_tbl$fdr < FDR_CUTOFF)),
  sprintf("Significant + overlap>=3 + FE>1: %d", nrow(sig_results)),
  "",
  "LITERATURE-ANCHORED VALIDATION",
  strrep("-", 60)
)

if (!is.null(rim101_check) || !is.null(sfp1_check) ||
    !is.null(gal4_check)   || !is.null(ppr1_check) ||
    !is.null(dal80_check)) {
  report_lines <- c(report_lines,
    "The following checks test whether motifs matching transcription",
    "factors with independently published Aspergillus/fungal regulatory",
    "roles show enrichment for the GO categories they are known to",
    "control, providing external validation of the RHIEPA pipeline",
    "independent of the GO annotations used to build it.",
    ""
  )
  if (!is.null(rim101_check)) {
    report_lines <- c(report_lines, sprintf(
      paste("MA0368.1 / RIM101 (PacC ortholog, alkaline pH response regulator,",
            "known to control cell wall remodeling and pathogenesis-related",
            "genes in Aspergillus) -> pathogenesis: overlap=%d/%d targets,",
            "%.2fx enrichment, FDR=%.3g"),
      rim101_check$overlap, rim101_check$target_set_size,
      rim101_check$fold_enrichment, rim101_check$fdr
    ), "")
  }
  if (!is.null(rim101_check2)) {
    report_lines <- c(report_lines, sprintf(
      paste("MA0368.1 / RIM101 -> cell wall/morphology: overlap=%d/%d targets,",
            "%.2fx enrichment, FDR=%.3g"),
      rim101_check2$overlap, rim101_check2$target_set_size,
      rim101_check2$fold_enrichment, rim101_check2$fdr
    ), "")
  }
  if (!is.null(sfp1_check)) {
    report_lines <- c(report_lines, sprintf(
      paste("MA0378.2 / SFP1 (ribosome biogenesis regulator, TORC1-responsive)",
            "-> ribosome biogenesis: overlap=%d/%d targets, %.2fx enrichment,",
            "FDR=%.3g"),
      sfp1_check$overlap, sfp1_check$target_set_size,
      sfp1_check$fold_enrichment, sfp1_check$fdr
    ), "")
  }
  if (!is.null(gal4_check)) {
    report_lines <- c(report_lines, sprintf(
      paste("MA0299.1 / GAL4 (galactose/carbon source regulator)",
            "-> carbohydrate metabolism: overlap=%d/%d targets,",
            "%.2fx enrichment, FDR=%.3g"),
      gal4_check$overlap, gal4_check$target_set_size,
      gal4_check$fold_enrichment, gal4_check$fdr
    ), "")
  }
  if (!is.null(ppr1_check)) {
    report_lines <- c(report_lines, sprintf(
      paste("MA2703.1 / PPR1 (pyrimidine pathway regulator)",
            "-> amino acid/vitamin metabolism: overlap=%d/%d targets,",
            "%.2fx enrichment, FDR=%.3g"),
      ppr1_check$overlap, ppr1_check$target_set_size,
      ppr1_check$fold_enrichment, ppr1_check$fdr
    ), "")
  }
  if (!is.null(dal80_check)) {
    report_lines <- c(report_lines, sprintf(
      paste("MA0289.1 / DAL80 (nitrogen catabolite repressor)",
            "-> nitrogen/transport: overlap=%d/%d targets,",
            "%.2fx enrichment, FDR=%.3g"),
      dal80_check$overlap, dal80_check$target_set_size,
      dal80_check$fold_enrichment, dal80_check$fdr
    ), "")
  }
}

report_lines <- c(report_lines,
  "OUTPUT FILES",
  strrep("-", 60),
  "  high_confidence_tf_motifs.tsv       - motifs passing IC/TomTom filter",
  "  tf_go_enrichment_full_results.tsv   - all pairwise test results",
  "  tf_go_enrichment_SIGNIFICANT.tsv    - FDR-significant, biologically",
  "                                        meaningful hits only",
  "  literature_anchored_checks.tsv      - specific TF/category checks",
  "  enrichment_volcano_plot.pdf/.png    - all tests, volcano style",
  "  top_significant_enrichments_barchart.pdf/.png - top 20 hits",
  ""
)

writeLines(report_lines, file.path(outdir, "SUMMARY_REPORT.txt"))
cat(paste(report_lines, collapse = "\n"))
cat(sprintf("\n\nAll outputs written to: %s\n", outdir))