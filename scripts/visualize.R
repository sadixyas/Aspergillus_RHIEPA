#!/usr/bin/env Rscript
# visualize.R
# -----------
# Produces:
#   1. TomTom summary plots
#   2. Motif information content plot
#   3. UMAP colored by Louvain clusters
#   4. UMAP colored by Aspergillus GO slim categories (one combined plot)
#   5. Individual UMAPs per GO slim category (highlighted vs grey)
#
# Uses the Aspergillus GO slim OBO from the cluster at:
#   /srv/projects/db/GO/current/goslim_aspergillus.obo
# And the full GO hierarchy from:
#   /srv/projects/db/GO/current/go-basic.obo
#
# Usage:
#   Rscript visualize.R motifs.meme tomtom.tsv rep.csv outdir \
#       og_map_filtered.tsv go.gaf.gz \
#       goslim_aspergillus.obo go-basic.obo

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript visualize.R motifs.meme tomtom.tsv rep.csv outdir",
       " [og_map] [gaf] [slim_obo] [go_basic_obo]\n", call. = FALSE)
}

meme_file  <- args[1]
tomtom_tsv <- args[2]
rep_csv    <- args[3]
outdir     <- args[4]
og_map     <- if (length(args) >= 5) args[5] else NULL
gaf_file   <- if (length(args) >= 6) args[6] else NULL
slim_obo   <- if (length(args) >= 7) args[7] else NULL
basic_obo  <- if (length(args) >= 8) args[8] else NULL

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

suppressMessages({
  library(readr); library(dplyr); library(stringr)
  library(tidyr); library(ggplot2)
})

# =============================================================
# FUNCTION: Parse an OBO file
# Returns list with:
#   $terms : tibble(id, name, namespace, is_slim, is_obsolete)
#   $parents: named list of character vectors (id -> parent ids)
# =============================================================
parse_obo <- function(path, slim_tag = NULL) {
  lines <- readLines(path, warn = FALSE)
  ids <- names <- nss <- character()
  slim <- obs <- logical()
  parents_list <- list()

  cur_id <- cur_name <- cur_ns <- NA_character_
  cur_slim <- cur_obs <- FALSE
  cur_parents <- character()
  in_term <- FALSE

  flush <- function() {
    if (!in_term || is.na(cur_id)) return()
    ids     <<- c(ids,     cur_id)
    names   <<- c(names,   cur_name)
    nss     <<- c(nss,     cur_ns)
    slim    <<- c(slim,    cur_slim)
    obs     <<- c(obs,     cur_obs)
    parents_list[[cur_id]] <<- cur_parents
  }

  for (line in lines) {
    if (line == "[Term]") {
      flush()
      cur_id <- cur_name <- cur_ns <- NA_character_
      cur_slim <- cur_obs <- FALSE
      cur_parents <- character()
      in_term <- TRUE
    } else if (in_term) {
      if      (startsWith(line, "id: "))
        cur_id <- sub("^id: ", "", line)
      else if (startsWith(line, "name: "))
        cur_name <- sub("^name: ", "", line)
      else if (startsWith(line, "namespace: "))
        cur_ns <- sub("^namespace: ", "", line)
      else if (line == "is_obsolete: true")
        cur_obs <- TRUE
      else if (!is.null(slim_tag) && startsWith(line, "subset:") &&
               grepl(slim_tag, line, fixed = TRUE))
        cur_slim <- TRUE
      else if (startsWith(line, "is_a: ")) {
        p <- str_extract(line, "GO:\\d+")
        if (!is.na(p)) cur_parents <- c(cur_parents, p)
      }
    }
  }
  flush()

  list(
    terms   = tibble(id=ids, name=names, namespace=nss,
                     is_slim=slim, is_obsolete=obs),
    parents = parents_list
  )
}

# =============================================================
# FUNCTION: Walk GO hierarchy to find nearest slim ancestor
# =============================================================
find_slim_ancestor <- function(go_id, parent_map, slim_set) {
  visited <- character()
  queue   <- go_id
  while (length(queue) > 0) {
    cur   <- queue[1]; queue <- queue[-1]
    if (cur %in% visited) next
    visited <- c(visited, cur)
    if (cur %in% slim_set) return(cur)
    pars <- parent_map[[cur]]
    if (length(pars) > 0) queue <- c(queue, pars)
  }
  NA_character_
}

# =============================================================
# 1. TomTom summary plots
# =============================================================
cat("=== 1. TomTom summary plots ===\n")
tt <- read_tsv(tomtom_tsv, comment="#", show_col_types=FALSE) %>%
  rename_with(~ gsub("-","_",.x))
stopifnot(all(c("Query_ID","Target_ID","E_value") %in% names(tt)))

best <- tt %>%
  group_by(Query_ID) %>%
  slice_min(E_value, n=1, with_ties=FALSE) %>% ungroup() %>%
  mutate(pass_strict = E_value < (1/256), pass_explore = E_value < 0.05)
n_strict <- sum(best$pass_strict)
cat(sprintf("  Motifs matched (E<1/256): %d / 256\n", n_strict))

ggsave(file.path(outdir,"tomtom_target_counts.pdf"),
  best %>% filter(pass_explore) %>% count(Target_ID,sort=TRUE) %>% head(20) %>%
  ggplot(aes(reorder(Target_ID,n),n)) + geom_col(fill="steelblue") +
  coord_flip() +
  labs(title=sprintf("Top TomTom targets (%d of 256 matched, E<1/256)",n_strict),
       x="Target motif", y="# matches") + theme_bw(),
  width=8, height=5)

ggsave(file.path(outdir,"tomtom_evalue_hist.pdf"),
  best %>% mutate(mlog10E=-log10(pmax(E_value,1e-15))) %>%
  ggplot(aes(mlog10E)) + geom_histogram(bins=40,fill="grey40") +
  geom_vline(xintercept=-log10(1/256),linetype="dashed",color="red") +
  geom_vline(xintercept=-log10(0.05), linetype="dotted",color="blue") +
  labs(title="TomTom E-value distribution",x="-log10(E-value)",y="count") +
  theme_bw(), width=6, height=4)
write_tsv(best, file.path(outdir,"tomtom_besthit_per_query.tsv"))

# =============================================================
# 2. Motif IC
# =============================================================
cat("=== 2. Motif information content ===\n")
col_ic    <- function(p) { p<-p[p>0]; sum(p*log2(p/0.25)) }
parse_meme <- function(path) {
  lines <- readLines(path,warn=FALSE); motifs <- list(); i <- 1
  while (i<=length(lines)) {
    if (startsWith(lines[i],"MOTIF ")) {
      nm <- str_trim(sub("^MOTIF\\s+","",lines[i])); j <- i+1
      while (j<=length(lines) && !str_detect(lines[j],"^letter-prob")) j<-j+1
      w  <- as.integer(str_match(lines[j],"w=\\s*([0-9]+)")[,2])
      mat <- do.call(rbind,lapply(lines[(j+1):(j+w)], function(x)
        as.numeric(strsplit(str_trim(x),"\\s+")[[1]])))
      mat <- mat/rowSums(mat); motifs[[nm]] <- mat; i <- j+w+1
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
  width=6, height=4)

# =============================================================
# 3. UMAP + Louvain
# =============================================================
cat("=== 3. UMAP + Louvain clustering ===\n")
rep      <- read_csv(rep_csv, show_col_types=FALSE)
num_cols <- rep %>% select(where(is.numeric)) %>% names()
cat(sprintf("  rep.csv: %d OGs x %d dimensions\n",nrow(rep),length(num_cols)))

if (length(num_cols) < 10) stop("rep.csv has <10 numeric columns")

X      <- rep %>% select(all_of(num_cols)) %>% as.matrix(); X[is.na(X)] <- 0
og_ids <- str_remove(rep$FILE,"\\.fa$")

for (pkg in c("uwot","igraph","FNN","RColorBrewer")) {
  if (!requireNamespace(pkg,quietly=TRUE))
    install.packages(pkg,repos="https://cran.r-project.org",quiet=TRUE)
}

cat("  Running UMAP...\n")
set.seed(42)
um    <- uwot::umap(X,n_neighbors=15,min_dist=0.2,metric="cosine",verbose=FALSE)
um_df <- tibble(OG=og_ids, UMAP1=um[,1], UMAP2=um[,2])

cat("  Running Louvain clustering...\n")
nn    <- FNN::get.knn(X,k=15)
edges <- do.call(rbind,lapply(seq_len(nrow(nn$nn.index)),
                               function(i) cbind(i,nn$nn.index[i,])))
g     <- igraph::simplify(igraph::graph_from_edgelist(edges,directed=FALSE))
louv  <- igraph::cluster_louvain(g)
um_df$cluster <- factor(igraph::membership(louv))
cat(sprintf("  Found %d Louvain clusters\n",length(unique(um_df$cluster))))

ggsave(file.path(outdir,"umap_louvain_clusters.pdf"),
  ggplot(um_df,aes(UMAP1,UMAP2,color=cluster)) +
  geom_point(alpha=0.6,size=0.8) +
  labs(title="UMAP of promoter representations",
       subtitle="Colors = unsupervised Louvain clusters",color="cluster") +
  theme_bw() +
  guides(color=guide_legend(override.aes=list(size=3,alpha=1),ncol=2)),
  width=9, height=6)
write_tsv(um_df, file.path(outdir,"umap_coordinates.tsv"))
cat("  Saved: umap_louvain_clusters.pdf\n")

# =============================================================
# 4. GO slim coloring
# =============================================================
if (!is.null(slim_obo) && !is.null(basic_obo) &&
    !is.null(og_map)   && !is.null(gaf_file)  &&
    file.exists(slim_obo) && file.exists(basic_obo) &&
    file.exists(og_map)   && file.exists(gaf_file)) {

  cat("=== 4. GO slim coloring (Aspergillus slim) ===\n")

  # Parse slim OBO to get the 37 BP slim term IDs and names
  cat("  Parsing Aspergillus GO slim OBO...\n")
  slim_obo_data <- parse_obo(slim_obo, slim_tag="goslim_aspergillus")
  slim_bp <- slim_obo_data$terms %>%
    filter(is_slim, !is_obsolete, namespace=="biological_process")
  cat(sprintf("  Slim BP terms: %d\n", nrow(slim_bp)))

  slim_ids   <- slim_bp$id
  slim_names <- setNames(slim_bp$name, slim_bp$id)

  # Parse go-basic.obo for hierarchy (parent map only — skip slim flag)
  cat("  Parsing go-basic.obo hierarchy...\n")
  basic_data <- parse_obo(basic_obo, slim_tag=NULL)
  parent_map <- basic_data$parents
  cat(sprintf("  GO terms with parent info: %d\n", length(parent_map)))

  # Parse GAF — biological process only
  cat("  Parsing GAF...\n")
  gaf_con   <- gzfile(gaf_file,"rt")
  gaf_lines <- readLines(gaf_con); close(gaf_con)
  gaf_lines <- gaf_lines[!startsWith(gaf_lines,"!")]
  gaf <- as.data.frame(do.call(rbind,strsplit(gaf_lines,"\t")),
                        stringsAsFactors=FALSE)
  gaf_bp <- gaf %>% select(gene_id=V2,go_term=V5,aspect=V9) %>%
    filter(aspect=="P") %>% distinct()
  cat(sprintf("  BP annotations: %d\n",nrow(gaf_bp)))

  # Map each unique GO term to its slim ancestor
  cat("  Mapping GO terms to slim ancestors...\n")
  unique_terms <- unique(gaf_bp$go_term)
  slim_map <- setNames(
    sapply(unique_terms, function(go_id) {
      # Check if it is itself a slim term
      if (go_id %in% slim_ids) return(go_id)
      find_slim_ancestor(go_id, parent_map, slim_ids)
    }),
    unique_terms
  )
  n_mapped <- sum(!is.na(slim_map))
  cat(sprintf("  GO terms mapped to slim ancestor: %d / %d\n",
              n_mapped, length(unique_terms)))

  # Priority ranking: specific biologically meaningful terms win
  # over broad catch-all terms. Lower number = higher priority.
  SLIM_PRIORITY <- c(
    "pathogenesis"                          = 1,
    "toxin metabolic process"               = 2,
    "secondary metabolic process"           = 3,
    "filamentous growth"                    = 4,
    "asexual sporulation"                   = 5,
    "sexual sporulation"                    = 6,
    "cell adhesion"                         = 7,
    "conjugation"                           = 8,
    "ribosome biogenesis"                   = 9,
    "translation"                           = 10,
    "cell cycle"                            = 11,
    "cytokinesis"                           = 12,
    "cellular respiration"                  = 13,
    "DNA metabolic process"                 = 14,
    "protein folding"                       = 15,
    "protein catabolic process"             = 16,
    "vitamin metabolic process"             = 17,
    "cellular amino acid metabolic process" = 18,
    "lipid metabolic process"               = 19,
    "carbohydrate metabolic process"        = 20,
    "vesicle-mediated transport"            = 21,
    "cellular homeostasis"                  = 22,
    "response to stress"                    = 23,
    "response to chemical"                  = 24,
    "signal transduction"                   = 25,
    "cytoskeleton organization"             = 26,
    "organelle organization"                = 27,
    "nucleus organization"                  = 28,
    "establishment of nucleus localization" = 29,
    "transport"                             = 30,
    "cellular protein modification process" = 31,
    "transcription, DNA-templated"          = 32,
    "RNA metabolic process"                 = 33,
    "developmental process"                 = 34,
    "regulation of biological process"      = 35
  )

  # Join slim ancestor to gene annotations, pick highest-priority slim term per gene
  gaf_slim <- gaf_bp %>%
    mutate(slim_id   = slim_map[go_term],
           slim_name = slim_names[slim_id]) %>%
    filter(!is.na(slim_name),
           slim_name != "biological_process") %>%
    mutate(priority = SLIM_PRIORITY[slim_name],
           priority = ifelse(is.na(priority), 99L, as.integer(priority))) %>%
    group_by(gene_id) %>%
    slice_min(priority, n=1, with_ties=FALSE) %>%
    ungroup() %>%
    select(gene_id, slim_name)
  cat(sprintf("  Genes with slim annotation: %d\n", nrow(gaf_slim)))

  # Load OG map, keep reference species, strip -T-p suffix
  og_map_df <- read_tsv(og_map, col_names=c("OG","species","gene_id"),
                         show_col_types=FALSE) %>%
    filter(str_detect(species,"FungiDB")) %>%
    mutate(gene_id_clean = str_remove(gene_id,"-T-p[0-9]+$"))

  # Join OG -> gene -> slim name
  og_slim <- og_map_df %>%
    left_join(gaf_slim, by=c("gene_id_clean"="gene_id")) %>%
    select(OG, slim_name) %>%
    mutate(slim_name = replace_na(slim_name, "unannotated"))

  um_go <- um_df %>%
    left_join(og_slim, by="OG") %>%
    mutate(slim_name = replace_na(slim_name, "unannotated"))

  n_ann <- sum(um_go$slim_name != "unannotated")
  cat(sprintf("  OGs with slim annotation: %d / %d\n", n_ann, nrow(um_go)))

  # Save full annotated coordinates
  write_tsv(um_go, file.path(outdir,"umap_coordinates_with_go.tsv"))

  # Category counts
  cat_counts <- um_go %>% filter(slim_name!="unannotated") %>%
    count(slim_name, sort=TRUE)
  write_tsv(cat_counts, file.path(outdir,"go_slim_category_counts.tsv"))
  cat("  GO slim categories found:\n"); print(cat_counts, n=40)

  # ── Combined UMAP colored by GO slim ────────────────────────
  cat_order  <- c(cat_counts$slim_name, "unannotated")
  um_go$slim_name <- factor(um_go$slim_name, levels=cat_order)
  n_cats <- length(cat_counts$slim_name)

  if (n_cats <= 12) {
    pal <- RColorBrewer::brewer.pal(max(3,n_cats),"Paired")[seq_len(n_cats)]
  } else {
    pal <- colorRampPalette(RColorBrewer::brewer.pal(12,"Paired"))(n_cats)
  }
  names(pal) <- cat_counts$slim_name

  um_bg   <- um_go %>% filter(slim_name=="unannotated")
  um_anno <- um_go %>% filter(slim_name!="unannotated")

  p_combined <- ggplot() +
    geom_point(data=um_bg,   aes(UMAP1,UMAP2),
               color="grey88", size=0.5, alpha=0.4) +
    geom_point(data=um_anno, aes(UMAP1,UMAP2,color=slim_name),
               size=1, alpha=0.85) +
    scale_color_manual(values=pal) +
    labs(title    = "UMAP of promoter representations",
         subtitle = sprintf("Aspergillus GO slim biological process (%d OGs annotated)",
                            n_ann),
         color="GO slim term") +
    theme_bw() +
    theme(legend.text=element_text(size=7)) +
    guides(color=guide_legend(override.aes=list(size=3,alpha=1),ncol=1))
  ggsave(file.path(outdir,"umap_go_slim_combined.pdf"),
         p_combined, width=12, height=7)
  cat("  Saved: umap_go_slim_combined.pdf\n")

  # ── Individual UMAP per GO slim category ────────────────────
  cat("  Generating individual UMAP per GO slim category...\n")
  indiv_dir <- file.path(outdir,"umap_per_go_slim")
  dir.create(indiv_dir, showWarnings=FALSE)

  # Base grey UMAP (all points)
  base_layer <- list(
    geom_point(data=um_go, aes(UMAP1,UMAP2),
               color="grey88", size=0.5, alpha=0.4)
  )

  for (cat_name in cat_counts$slim_name) {
    hi <- um_go %>% filter(slim_name==cat_name)
    n_hi <- nrow(hi)
    cat_color <- pal[cat_name]

    p <- ggplot() +
      base_layer +
      geom_point(data=hi, aes(UMAP1,UMAP2),
                 color=cat_color, size=1.2, alpha=0.9) +
      labs(title    = cat_name,
           subtitle = sprintf("%d OGs", n_hi),
           x="UMAP1", y="UMAP2") +
      theme_bw() +
      theme(plot.title=element_text(size=12, face="bold"),
            plot.subtitle=element_text(size=9))

    # Sanitize filename
    safe_name <- str_replace_all(cat_name, "[^A-Za-z0-9_]", "_")
    ggsave(file.path(indiv_dir, paste0(safe_name,".pdf")),
           p, width=6, height=5)
  }
  cat(sprintf("  Saved %d individual UMAPs to: %s\n",
              nrow(cat_counts), indiv_dir))

  # ── Multi-panel figure (all categories in one PDF) ──────────
  cat("  Generating multi-panel figure...\n")
  if (!requireNamespace("patchwork", quietly=TRUE))
    install.packages("patchwork", repos="https://cran.r-project.org", quiet=TRUE)
  library(patchwork)

  plot_list <- lapply(cat_counts$slim_name, function(cat_name) {
    hi <- um_go %>% filter(slim_name==cat_name)
    ggplot() +
      geom_point(data=um_go, aes(UMAP1,UMAP2),
                 color="grey88", size=0.3, alpha=0.3) +
      geom_point(data=hi, aes(UMAP1,UMAP2),
                 color=pal[cat_name], size=0.8, alpha=0.9) +
      labs(title=cat_name,
           subtitle=sprintf("n=%d", nrow(hi))) +
      theme_bw(base_size=8) +
      theme(plot.title=element_text(size=7,face="bold"),
            plot.subtitle=element_text(size=6),
            axis.title=element_blank(),
            axis.text=element_text(size=5))
  })

  ncols <- 5
  nrows <- ceiling(length(plot_list)/ncols)
  combined <- wrap_plots(plot_list, ncol=ncols)
  ggsave(file.path(outdir,"umap_go_slim_multipanel.pdf"),
         combined, width=ncols*4, height=nrows*3.5)
  cat("  Saved: umap_go_slim_multipanel.pdf\n")

} else {
  cat("  Skipping GO slim coloring: one or more input files not found\n")
  cat(sprintf("  slim_obo: %s\n", ifelse(is.null(slim_obo),"NULL",slim_obo)))
  cat(sprintf("  basic_obo: %s\n", ifelse(is.null(basic_obo),"NULL",basic_obo)))
  cat(sprintf("  og_map: %s\n", ifelse(is.null(og_map),"NULL",og_map)))
  cat(sprintf("  gaf_file: %s\n", ifelse(is.null(gaf_file),"NULL",gaf_file)))
}

cat("\nDONE - results in:", outdir, "\n")