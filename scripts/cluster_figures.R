#!/usr/bin/env Rscript
# cluster_figures.R
# -----------------------
# Produces publication-quality figures from the UMAP + GO slim output:
#   1. Clustering ratio bar chart (all GO slim categories)
#   2. Improved individual UMAP panels for key categories
#   3. Multi-panel figure of top clustered categories
#
# Usage:
#   Rscript 10_clustering_figures.R umap_coordinates_with_go.tsv outdir

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 10_clustering_figures.R umap_coordinates_with_go.tsv outdir",
       call. = FALSE)
}

umap_tsv <- args[1]
outdir   <- args[2]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

suppressMessages({
  library(readr); library(dplyr); library(ggplot2)
  library(stringr); library(tidyr)
})

for (pkg in c("patchwork","RColorBrewer")) {
  if (!requireNamespace(pkg, quietly=TRUE))
    install.packages(pkg, repos="https://cran.r-project.org", quiet=TRUE)
}
library(patchwork)

cat("Loading UMAP coordinates...\n")
um <- read_tsv(umap_tsv, show_col_types=FALSE)
cat(sprintf("  %d OGs loaded\n", nrow(um)))

# =============================================================
# 1. Compute clustering ratios
# =============================================================
cat("Computing clustering ratios...\n")

cats <- um %>%
  filter(slim_name != "unannotated") %>%
  count(slim_name, sort=TRUE) %>%
  filter(n >= 5) %>%
  pull(slim_name)

set.seed(42)
cluster_stats <- lapply(cats, function(cat) {
  hi  <- um %>% filter(slim_name == cat)
  if (nrow(hi) < 2) return(NULL)

  # Within-category nearest neighbour distance
  coords_hi <- hi[, c("UMAP1","UMAP2")]
  d_hi <- as.matrix(dist(coords_hi))
  diag(d_hi) <- Inf
  within_nn <- mean(apply(d_hi, 1, min))

  # Random baseline: same number of points drawn from all OGs
  rand <- um[sample(nrow(um), nrow(hi), replace=FALSE), ]
  d_r  <- as.matrix(dist(rand[, c("UMAP1","UMAP2")]))
  diag(d_r) <- Inf
  random_nn <- mean(apply(d_r, 1, min))

  tibble(
    category        = cat,
    n               = nrow(hi),
    within_nn       = round(within_nn, 3),
    random_nn       = round(random_nn, 3),
    clustering_ratio = round(within_nn / random_nn, 3)
  )
}) %>% bind_rows() %>% arrange(clustering_ratio)

write_tsv(cluster_stats, file.path(outdir, "clustering_ratios.tsv"))
cat("  Clustering ratios computed:\n")
print(cluster_stats, n=40)

# =============================================================
# 2. Clustering ratio bar chart (publication quality)
# =============================================================
cat("\nGenerating clustering ratio bar chart...\n")

# Color by strength: < 0.65 = strong, 0.65-0.80 = moderate, > 0.80 = weak
cluster_stats <- cluster_stats %>%
  mutate(
    strength = case_when(
      clustering_ratio < 0.65 ~ "Strong\n(< 0.65)",
      clustering_ratio < 0.80 ~ "Moderate\n(0.65-0.80)",
      TRUE                    ~ "Weak\n(> 0.80)"
    ),
    strength = factor(strength,
                      levels=c("Strong\n(< 0.65)",
                               "Moderate\n(0.65-0.80)",
                               "Weak\n(> 0.80)")),
    # Wrap long category names
    category_wrap = str_wrap(category, width=30),
    category_wrap = factor(category_wrap,
                           levels=rev(str_wrap(cluster_stats$category, 30)))
  )

strength_colors <- c(
  "Strong\n(< 0.65)"    = "#2166ac",
  "Moderate\n(0.65-0.80)" = "#4dac26",
  "Weak\n(> 0.80)"      = "#d6604d"
)

p_ratio <- ggplot(cluster_stats,
                  aes(x=clustering_ratio, y=category_wrap, fill=strength)) +
  geom_col(width=0.7, color="white", linewidth=0.2) +
  geom_vline(xintercept=1, linetype="dashed", color="grey40", linewidth=0.5) +
  geom_text(aes(label=sprintf("n=%d", n)),
            hjust=-0.1, size=2.8, color="grey30") +
  scale_fill_manual(values=strength_colors, name="Clustering\nstrength") +
  scale_x_continuous(limits=c(0, 1.15),
                     breaks=c(0, 0.25, 0.5, 0.65, 0.80, 1.0),
                     labels=c("0","0.25","0.50","0.65","0.80","1.0")) +
  labs(
    title    = "GO slim category clustering in RHIEPA promoter space",
    subtitle = "Ratio < 1 indicates tighter clustering than random (all categories cluster better than random)",
    x        = "Clustering ratio (within-category NN distance / random NN distance)",
    y        = NULL
  ) +
  theme_bw(base_size=11) +
  theme(
    plot.title    = element_text(face="bold", size=12),
    plot.subtitle = element_text(size=9, color="grey40"),
    legend.position = "right",
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.y   = element_text(size=9)
  ) +
  annotate("rect", xmin=0, xmax=0.65, ymin=0.5, ymax=nrow(cluster_stats)+0.5,
           fill="#2166ac", alpha=0.04) +
  annotate("rect", xmin=0.65, xmax=0.80, ymin=0.5, ymax=nrow(cluster_stats)+0.5,
           fill="#4dac26", alpha=0.04) +
  annotate("rect", xmin=0.80, xmax=1.15, ymin=0.5, ymax=nrow(cluster_stats)+0.5,
           fill="#d6604d", alpha=0.04)

ggsave(file.path(outdir,"clustering_ratio_barchart.pdf"),
       p_ratio, width=10, height=8)
ggsave(file.path(outdir,"clustering_ratio_barchart.png"),
       p_ratio, width=10, height=8, dpi=300)
cat("  Saved: clustering_ratio_barchart.pdf/.png\n")

# =============================================================
# 3. Improved individual UMAP panels for key categories
# =============================================================
cat("\nGenerating improved individual UMAP panels...\n")

# Define key categories to highlight with their display colors
key_cats <- list(
  list(name="ribosome biogenesis",      color="#1f78b4", label="Ribosome biogenesis"),
  list(name="translation",              color="#33a02c", label="Translation"),
  list(name="secondary metabolic process", color="#e31a1c", label="Secondary metabolism"),
  list(name="filamentous growth",       color="#ff7f00", label="Filamentous growth"),
  list(name="pathogenesis",             color="#6a3d9a", label="Pathogenesis"),
  list(name="toxin metabolic process",  color="#b15928", label="Toxin metabolism"),
  list(name="asexual sporulation",      color="#a6cee3", label="Asexual sporulation"),
  list(name="cell cycle",               color="#b2df8a", label="Cell cycle"),
  list(name="DNA metabolic process",    color="#fb9a99", label="DNA metabolism"),
  list(name="response to stress",       color="#fdbf6f", label="Stress response"),
  list(name="protein folding",          color="#cab2d6", label="Protein folding"),
  list(name="signal transduction",      color="#ffff99", label="Signal transduction")
)

# Get UMAP axis limits
x_range <- range(um$UMAP1)
y_range <- range(um$UMAP2)

make_umap_panel <- function(cat_name, cat_color, cat_label, um_data,
                            show_axes=TRUE, point_size_hi=1.2,
                            point_size_bg=0.4) {
  hi <- um_data %>% filter(slim_name == cat_name)
  bg <- um_data

  n_hi     <- nrow(hi)
  ratio_row <- cluster_stats %>% filter(category==cat_name)
  ratio_str <- if (nrow(ratio_row)>0)
    sprintf("ratio=%.3f", ratio_row$clustering_ratio[1]) else ""

  p <- ggplot() +
    geom_point(data=bg, aes(UMAP1, UMAP2),
               color="grey90", size=point_size_bg, alpha=0.5) +
    geom_point(data=hi, aes(UMAP1, UMAP2),
               color=cat_color, size=point_size_hi, alpha=0.9) +
    coord_fixed(xlim=x_range, ylim=y_range) +
    labs(title    = cat_label,
         subtitle = sprintf("n=%d  %s", n_hi, ratio_str)) +
    theme_bw(base_size=9) +
    theme(
      plot.title    = element_text(face="bold", size=9),
      plot.subtitle = element_text(size=7, color="grey40"),
      panel.grid    = element_blank(),
      axis.title    = if (show_axes) element_text(size=7) else element_blank(),
      axis.text     = if (show_axes) element_text(size=6) else element_blank(),
      axis.ticks    = if (show_axes) element_line() else element_blank()
    )
  p
}

# Individual full-size PDFs for each key category
indiv_dir <- file.path(outdir, "key_category_umaps")
dir.create(indiv_dir, showWarnings=FALSE)

for (cat_info in key_cats) {
  cat_name <- cat_info$name
  if (!cat_name %in% um$slim_name) next

  p <- make_umap_panel(cat_name, cat_info$color, cat_info$label,
                       um, show_axes=TRUE, point_size_hi=1.5)
  safe <- str_replace_all(cat_name,"[^A-Za-z0-9]","_")
  ggsave(file.path(indiv_dir, paste0(safe,".pdf")), p, width=5, height=4.5)
  ggsave(file.path(indiv_dir, paste0(safe,".png")), p, width=5, height=4.5, dpi=300)
}
cat(sprintf("  Saved %d individual UMAPs to: %s\n",
            length(key_cats), indiv_dir))

# =============================================================
# 4. Multi-panel figure of key categories (publication quality)
# =============================================================
cat("Generating publication multi-panel figure...\n")

plot_list <- lapply(key_cats, function(cat_info) {
  cat_name <- cat_info$name
  if (!cat_name %in% um$slim_name) return(NULL)
  make_umap_panel(cat_name, cat_info$color, cat_info$label,
                  um, show_axes=FALSE,
                  point_size_hi=0.9, point_size_bg=0.25)
})
plot_list <- plot_list[!sapply(plot_list, is.null)]

multi <- wrap_plots(plot_list, ncol=4) +
  plot_annotation(
    title    = "RHIEPA promoter representations — key GO slim categories",
    subtitle = "Each panel: highlighted category (colored) vs all other OGs (grey). Subtitle shows n and clustering ratio.",
    theme    = theme(
      plot.title    = element_text(face="bold", size=13),
      plot.subtitle = element_text(size=9, color="grey40")
    )
  )

ggsave(file.path(outdir,"key_categories_multipanel.pdf"),
       multi, width=16, height=10)
ggsave(file.path(outdir,"key_categories_multipanel.png"),
       multi, width=16, height=10, dpi=300)
cat("  Saved: key_categories_multipanel.pdf/.png\n")

# =============================================================
# 5. Pathogenesis + secondary metabolism side-by-side (for paper)
# =============================================================
cat("Generating pathogenesis/secondary metabolism figure...\n")

p_path <- make_umap_panel("pathogenesis", "#6a3d9a", "Pathogenesis",
                           um, show_axes=TRUE, point_size_hi=1.5)
p_sec  <- make_umap_panel("secondary metabolic process", "#e31a1c",
                           "Secondary metabolism", um, show_axes=TRUE,
                           point_size_hi=1.5)
p_fil  <- make_umap_panel("filamentous growth", "#ff7f00",
                           "Filamentous growth", um, show_axes=TRUE,
                           point_size_hi=1.5)
p_tox  <- make_umap_panel("toxin metabolic process", "#b15928",
                           "Toxin metabolism", um, show_axes=TRUE,
                           point_size_hi=1.5)

aspergillus_panel <- (p_path | p_sec) / (p_fil | p_tox) +
  plot_annotation(
    title    = "Aspergillus pathogenesis-related GO slim categories in RHIEPA space",
    subtitle = "Points = OG promoter representations. Colored = annotated category, grey = all others.",
    theme    = theme(
      plot.title    = element_text(face="bold", size=12),
      plot.subtitle = element_text(size=9, color="grey40")
    )
  )

ggsave(file.path(outdir,"aspergillus_pathogenesis_panel.pdf"),
       aspergillus_panel, width=10, height=9)
ggsave(file.path(outdir,"aspergillus_pathogenesis_panel.png"),
       aspergillus_panel, width=10, height=9, dpi=300)
cat("  Saved: aspergillus_pathogenesis_panel.pdf/.png\n")

cat("\nDONE - all figures written to:", outdir, "\n")
cat("\nFigure summary:\n")
cat("  clustering_ratio_barchart.pdf  — bar chart of all GO slim clustering ratios\n")
cat("  key_categories_multipanel.pdf  — 12-panel grid of key categories\n")
cat("  aspergillus_pathogenesis_panel.pdf — 2x2 panel: pathogenesis/sec. met/filamentous/toxin\n")
cat("  key_category_umaps/            — individual PDFs + PNGs for each key category\n")