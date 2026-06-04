#!/usr/bin/env Rscript
# visualize.R
# -----------
# Produces:
#   1. TomTom summary plots
#   2. Motif information content plots
#   3. UMAP colored by Louvain clusters (unsupervised)
#   4. UMAP colored by broad GO biological process category
#      (categories assigned by keyword matching on GO term names
#       from go-basic.obo — no GO slim file needed)
#
# Usage:
#   Rscript visualize.R motifs.meme tomtom.tsv rep.csv outdir \
#       og_map_filtered.tsv go.gaf.gz go-basic.obo

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript visualize.R motifs.meme tomtom.tsv rep.csv outdir",
       " [og_map.tsv] [go.gaf.gz] [go-basic.obo]\n", call. = FALSE)
}

meme_file  <- args[1]
tomtom_tsv <- args[2]
rep_csv    <- args[3]
outdir     <- args[4]
og_map     <- if (length(args) >= 5) args[5] else NULL
gaf_file   <- if (length(args) >= 6) args[6] else NULL
obo_file   <- if (length(args) >= 7) args[7] else NULL

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

suppressMessages({
  library(readr); library(dplyr); library(stringr)
  library(tidyr); library(ggplot2)
})

# =============================================================
# FUNCTION: Parse go-basic.obo -> tibble of id, name, namespace
# =============================================================
parse_obo_names <- function(obo_path) {
  cat("  Parsing GO term names from OBO file...\n")
  lines <- readLines(obo_path, warn = FALSE)

  ids   <- character()
  names <- character()
  nss   <- character()
  obs   <- logical()

  current_id  <- NA
  current_name <- NA
  current_ns  <- NA
  current_obs <- FALSE
  in_term     <- FALSE

  for (line in lines) {
    if (line == "[Term]") {
      if (in_term && !is.na(current_id)) {
        ids   <- c(ids,   current_id)
        names <- c(names, current_name)
        nss   <- c(nss,   current_ns)
        obs   <- c(obs,   current_obs)
      }
      current_id   <- NA
      current_name <- NA
      current_ns   <- NA
      current_obs  <- FALSE
      in_term      <- TRUE
    } else if (in_term) {
      if (startsWith(line, "id: "))
        current_id <- sub("^id: ", "", line)
      else if (startsWith(line, "name: "))
        current_name <- sub("^name: ", "", line)
      else if (startsWith(line, "namespace: "))
        current_ns <- sub("^namespace: ", "", line)
      else if (line == "is_obsolete: true")
        current_obs <- TRUE
    }
  }
  # last term
  if (in_term && !is.na(current_id)) {
    ids   <- c(ids,   current_id)
    names <- c(names, current_name)
    nss   <- c(nss,   current_ns)
    obs   <- c(obs,   current_obs)
  }

  tibble(id = ids, name = names, namespace = nss, is_obsolete = obs) %>%
    filter(!is_obsolete, namespace == "biological_process")
}

# =============================================================
# FUNCTION: Assign broad category from GO term name by keywords
# Rules are checked in order — first match wins
# =============================================================
assign_category <- function(term_name) {
  n <- tolower(term_name)

  if (str_detect(n, "ribosom|translation(?! regul)|translational"))
    return("Ribosome/Translation")
  if (str_detect(n, "cell cycle|mitosis|meiosis|cytokinesis|chromosome segreg"))
    return("Cell cycle")
  if (str_detect(n, "dna replicat|dna repair|dna damage|dna integr"))
    return("DNA replication/repair")
  if (str_detect(n, "transcription|mrna|rna polymerase|gene expression"))
    return("Transcription/RNA")
  if (str_detect(n, "splicing|rna process|rna metabol"))
    return("RNA processing")
  if (str_detect(n, "transport|transmembrane|import|export|vesicle|secreti"))
    return("Transport/Secretion")
  if (str_detect(n, "amino acid|protein folding|proteolysis|ubiquitin|proteasom"))
    return("Protein metabolism")
  if (str_detect(n, "carbohydrate|glycolysis|gluconeogenesis|sugar|polysaccharide"))
    return("Carbohydrate metabolism")
  if (str_detect(n, "lipid|fatty acid|sterol|ergosterol"))
    return("Lipid metabolism")
  if (str_detect(n, "secondary metabol|polyketide|terpene|alkaloid|mycotoxin|aflatoxin|fumonisin"))
    return("Secondary metabolism")
  if (str_detect(n, "pathogen|virulence|infect|host|effector"))
    return("Pathogenesis")
  if (str_detect(n, "stress|heat shock|oxidative|osmotic|response to"))
    return("Stress response")
  if (str_detect(n, "mitochondri|respiration|oxidative phosphoryl|electron transport"))
    return("Mitochondria/Respiration")
  if (str_detect(n, "cell wall|chitin|glucan|hyphal|conidiat|sporulat|spore"))
    return("Cell wall/Morphology")
  if (str_detect(n, "mating|pheromone|sexual|asexual|reproduction|fertiliz"))
    return("Reproduction")
  if (str_detect(n, "nitrogen|sulfur|phosphate|ion|metal|iron|zinc|copper"))
    return("Nutrient/Ion metabolism")
  if (str_detect(n, "signal|kinase|phosphorylation|mapk|camp|gtpase"))
    return("Signaling")
  if (str_detect(n, "autophagy|vacuole|lysosome|endosome|golgi|er "))
    return("Membrane/Organelle")
  if (str_detect(n, "metabol|biosynth|catabolic|biosynthetic"))
    return("Other metabolism")

  return("Other/Unknown")
}

# =============================================================
# 1. TomTom summary plots
# =============================================================
cat("=== 1. TomTom summary plots ===\n")
tt <- read_tsv(tomtom_tsv, comment = "#", show_col_types = FALSE) %>%
  rename_with(~ gsub("-", "_", .x))
stopifnot(all(c("Query_ID", "Target_ID", "E_value") %in% names(tt)))

best <- tt %>%
  group_by(Query_ID) %>%
  slice_min(E_value, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(pass_strict  = E_value < (1/256),
         pass_explore = E_value < 0.05)

n_strict <- sum(best$pass_strict)
cat(sprintf("  Motifs matched (E < 1/256): %d / 256\n", n_strict))

top_targets <- best %>% filter(pass_explore) %>%
  count(Target_ID, sort = TRUE) %>% head(20)

p_counts <- ggplot(top_targets, aes(reorder(Target_ID, n), n)) +
  geom_col(fill = "steelblue") + coord_flip() +
  labs(title   = sprintf("Top TomTom targets: %d of 256 motifs matched (E<1/256)",
                          n_strict),
       x = "Target motif", y = "# significant motif matches") +
  theme_bw()
ggsave(file.path(outdir, "tomtom_target_counts.pdf"), p_counts, w=8, h=5)

p_e <- best %>% mutate(mlog10E = -log10(pmax(E_value, 1e-15))) %>%
  ggplot(aes(mlog10E)) +
  geom_histogram(bins=40, fill="grey40") +
  geom_vline(xintercept=-log10(1/256), linetype="dashed", color="red") +
  geom_vline(xintercept=-log10(0.05),  linetype="dotted", color="blue") +
  labs(title="TomTom significance", subtitle="Red=E<1/256 | Blue=E<0.05",
       x="-log10(E-value)", y="count") + theme_bw()
ggsave(file.path(outdir, "tomtom_evalue_hist.pdf"), p_e, w=6, h=4)
write_tsv(best, file.path(outdir, "tomtom_besthit_per_query.tsv"))

# =============================================================
# 2. Motif information content
# =============================================================
cat("=== 2. Motif information content ===\n")
col_ic <- function(p) { p <- p[p>0]; sum(p * log2(p/0.25)) }

parse_meme <- function(path) {
  lines <- readLines(path, warn=FALSE); motifs <- list(); i <- 1
  while (i <= length(lines)) {
    if (startsWith(lines[i], "MOTIF ")) {
      name <- str_trim(sub("^MOTIF\\s+","",lines[i]))
      j <- i+1
      while (j<=length(lines) && !str_detect(lines[j],"^letter-probability")) j<-j+1
      w <- as.integer(str_match(lines[j],"w=\\s*([0-9]+)")[,2])
      mat <- do.call(rbind, lapply(lines[(j+1):(j+w)], function(x)
        as.numeric(strsplit(str_trim(x),"\\s+")[[1]])))
      mat <- mat/rowSums(mat); motifs[[name]] <- mat; i <- j+w+1
    } else i <- i+1
  }
  motifs
}

pwms <- parse_meme(meme_file)
ic_tbl <- tibble(motif=names(pwms),
  ic_total=sapply(pwms,function(m) sum(apply(m,1,col_ic))),
  ic_mean =sapply(pwms,function(m) mean(apply(m,1,col_ic)))) %>%
  arrange(desc(ic_total))
write_tsv(ic_tbl, file.path(outdir,"motif_information_content.tsv"))
ggsave(file.path(outdir,"motif_ic_hist.pdf"),
  ggplot(ic_tbl,aes(ic_total))+geom_histogram(bins=40,fill="grey40")+
  labs(title="Motif IC distribution",x="IC total (bits)",y="count")+theme_bw(),
  w=6,h=4)

# =============================================================
# 3. UMAP
# =============================================================
cat("=== 3. UMAP ===\n")
rep      <- read_csv(rep_csv, show_col_types=FALSE)
num_cols <- rep %>% select(where(is.numeric)) %>% names()
cat(sprintf("  rep.csv: %d OGs x %d dimensions\n", nrow(rep), length(num_cols)))

if (length(num_cols) < 10) {
  warning("rep.csv has <10 numeric columns; skipping UMAP.")
} else {
  X      <- rep %>% select(all_of(num_cols)) %>% as.matrix()
  X[is.na(X)] <- 0
  og_ids <- str_remove(rep$FILE, "\\.fa$")

  for (pkg in c("uwot","igraph","FNN","RColorBrewer")) {
    if (!requireNamespace(pkg, quietly=TRUE))
      install.packages(pkg, repos="https://cran.r-project.org", quiet=TRUE)
  }

  cat("  Running UMAP...\n")
  set.seed(42)
  um    <- uwot::umap(X, n_neighbors=15, min_dist=0.2,
                      metric="cosine", verbose=FALSE)
  um_df <- tibble(OG=og_ids, UMAP1=um[,1], UMAP2=um[,2])

  # Louvain clustering
  cat("  Running Louvain clustering...\n")
  nn    <- FNN::get.knn(X, k=15)
  edges <- do.call(rbind, lapply(seq_len(nrow(nn$nn.index)),
                                  function(i) cbind(i, nn$nn.index[i,])))
  g       <- igraph::simplify(igraph::graph_from_edgelist(edges, directed=FALSE))
  louv    <- igraph::cluster_louvain(g)
  um_df$cluster <- factor(igraph::membership(louv))
  cat(sprintf("  Found %d Louvain clusters\n", length(unique(um_df$cluster))))

  p_louv <- ggplot(um_df, aes(UMAP1, UMAP2, color=cluster)) +
    geom_point(alpha=0.6, size=0.8) +
    labs(title="UMAP of promoter representations",
         subtitle="Colors = unsupervised Louvain clusters",
         color="cluster") + theme_bw() +
    guides(color=guide_legend(override.aes=list(size=3,alpha=1), ncol=2))
  ggsave(file.path(outdir,"umap_louvain_clusters.pdf"), p_louv, w=9, h=6)
  cat("  Saved: umap_louvain_clusters.pdf\n")
  write_tsv(um_df, file.path(outdir,"umap_coordinates.tsv"))

  # =============================================================
  # 4. GO keyword coloring
  # =============================================================
  if (!is.null(og_map) && !is.null(gaf_file) && !is.null(obo_file) &&
      file.exists(og_map) && file.exists(gaf_file) && file.exists(obo_file)) {

    cat("=== 4. GO category coloring (keyword matching) ===\n")

    # Parse GO term names from OBO
    go_terms <- parse_obo_names(obo_file)
    cat(sprintf("  Loaded %d biological process GO terms\n", nrow(go_terms)))

    # Assign broad category to each GO term by keyword matching
    go_terms <- go_terms %>%
      mutate(category = sapply(name, assign_category))

    cat("  Category distribution in GO terms:\n")
    print(go_terms %>% count(category, sort=TRUE) %>% head(15))

    # Parse GAF (biological process only)
    cat("  Parsing GAF...\n")
    gaf_con   <- gzfile(gaf_file, "rt")
    gaf_lines <- readLines(gaf_con); close(gaf_con)
    gaf_lines <- gaf_lines[!startsWith(gaf_lines, "!")]
    gaf       <- as.data.frame(do.call(rbind, strsplit(gaf_lines, "\t")),
                                stringsAsFactors=FALSE)
    gaf_bp <- gaf %>%
      select(gene_id=V2, go_term=V5, aspect=V9) %>%
      filter(aspect=="P") %>%
      distinct()
    cat(sprintf("  Biological process annotations: %d\n", nrow(gaf_bp)))

    # Join GO term names -> categories
    gaf_cat <- gaf_bp %>%
      left_join(go_terms %>% select(id, category), by=c("go_term"="id")) %>%
      filter(!is.na(category), category != "Other/Unknown") %>%
      group_by(gene_id) %>%
      # pick most specific (non-Other) category per gene
      count(gene_id, category, sort=TRUE) %>%
      slice(1) %>% ungroup() %>%
      select(gene_id, category)

    cat(sprintf("  Genes with non-Other GO category: %d\n", nrow(gaf_cat)))

    # Load OG map, keep reference, strip -T-p1 suffix
    og_map_df <- read_tsv(og_map, col_names=c("OG","species","gene_id"),
                           show_col_types=FALSE) %>%
      filter(str_detect(species, "FungiDB")) %>%
      mutate(gene_id_clean = str_remove(gene_id, "-T-p[0-9]+$"))

    # Join OG -> gene -> category
    og_go <- og_map_df %>%
      left_join(gaf_cat, by=c("gene_id_clean"="gene_id")) %>%
      select(OG, category) %>%
      mutate(category = replace_na(category, "Other/Unknown"))

    um_go <- um_df %>%
      left_join(og_go, by="OG") %>%
      mutate(category = replace_na(category, "Other/Unknown"))

    n_ann <- sum(um_go$category != "Other/Unknown")
    cat(sprintf("  OGs with GO category: %d / %d\n", n_ann, nrow(um_go)))

    # Write category counts
    cat_counts <- um_go %>% filter(category!="Other/Unknown") %>%
      count(category, sort=TRUE)
    write_tsv(cat_counts, file.path(outdir,"go_category_counts.tsv"))
    cat("  Categories found:\n"); print(cat_counts)

    # Order by frequency for legend
    cat_order  <- c(cat_counts$category, "Other/Unknown")
    um_go$category <- factor(um_go$category, levels=cat_order)

    um_bg   <- um_go %>% filter(category=="Other/Unknown")
    um_anno <- um_go %>% filter(category!="Other/Unknown")
    n_cats  <- length(unique(um_anno$category))

    # Color palette
    if (n_cats <= 12) {
      pal <- RColorBrewer::brewer.pal(max(3,n_cats), "Paired")[1:n_cats]
    } else {
      pal <- colorRampPalette(RColorBrewer::brewer.pal(12,"Paired"))(n_cats)
    }

    p_go <- ggplot() +
      geom_point(data=um_bg,   aes(UMAP1,UMAP2),
                 color="grey88", size=0.5, alpha=0.4) +
      geom_point(data=um_anno, aes(UMAP1,UMAP2, color=category),
                 size=1, alpha=0.85) +
      scale_color_manual(values=pal) +
      labs(title    = "UMAP of promoter representations",
           subtitle = sprintf("Colors = broad GO biological process (%d OGs annotated)",
                              n_ann),
           color = "GO category") +
      theme_bw() +
      theme(legend.text=element_text(size=8)) +
      guides(color=guide_legend(override.aes=list(size=3,alpha=1), ncol=1))

    ggsave(file.path(outdir,"umap_go_categories.pdf"), p_go, w=12, h=7)
    write_tsv(um_go, file.path(outdir,"umap_coordinates_with_go.tsv"))
    cat("  Saved: umap_go_categories.pdf\n")

  } else {
    cat("  Skipping GO coloring: og_map, gaf_file, or obo_file not provided\n")
  }
}

cat("\nDONE - results in:", outdir, "\n")