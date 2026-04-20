#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript parse_tomtom.R tomtom.tsv outdir", call.=FALSE)
}
infile <- args[1]
outdir <- args[2]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

suppressWarnings({
  suppressMessages({
    library(readr)
    library(dplyr)
    library(stringr)
  })
})

# TomTom TSV has comment lines starting with '#'
tt <- read_tsv(infile, comment = "#", show_col_types = FALSE)

# Common columns include: Query_ID, Target_ID, Optimal_offset, p-value, E-value, q-value, Overlap, Query_consensus, Target_consensus, Orientation
# Some builds use slightly different names; standardize a bit:
names(tt) <- str_replace_all(names(tt), "[ -]", "_")

# Pick best hit per query by q-value then E-value then p-value
best <- tt %>%
  filter(!is.na(Query_ID), !is.na(Target_ID)) %>%
  arrange(Query_ID, q_value, E_value, p_value) %>%
  group_by(Query_ID) %>%
  slice(1) %>%
  ungroup()

# Split Target_ID like "MA0271.1 ARG80" when present
best <- best %>%
  mutate(
    target_acc = str_extract(Target_ID, "^MA\\d+\\.\\d+"),
    target_name = str_trim(str_remove(Target_ID, "^MA\\d+\\.\\d+\\s*"))
  )

out <- file.path(outdir, "motif_annotations.tsv")
write_tsv(best, out)

cat("Wrote:", out, "\n")
cat("Queries with a hit:", nrow(best), "\n")
cat("Example:\n")
print(head(best %>% select(Query_ID, Target_ID, target_acc, target_name, q_value, E_value, p_value), 10))