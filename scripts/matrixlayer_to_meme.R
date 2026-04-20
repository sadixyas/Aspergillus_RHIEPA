#!/usr/bin/env Rscript
# ============================================================
# 00_matrixlayer_to_meme_ACGTlines.R
# Convert Moses matrix_layer.txt (A/C/G/T lines with 15 numbers each)
# into MEME motif format.
#
# Expected per motif block:
#   >matrix0 length=15  thr=... sca=... pool=... interactions: ...
#   A <15 nums>
#   C <15 nums>
#   G <15 nums>
#   T <15 nums>
#
# Usage:
#   Rscript matrixlayer_to_meme.R matrix_layer.txt motifs.meme [filter_len]
#
# Example:
#   Rscript matrixlayer_to_meme.R promoter_rh_*_matrix_layer.txt motifs.meme 15
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript matrixlayer_to_meme.R matrix_layer.txt motifs.meme [filter_len]", call. = FALSE)
}

in_path   <- args[1]
out_path  <- args[2]
filter_len <- if (length(args) >= 3) as.integer(args[3]) else 15L

if (!file.exists(in_path)) stop("Input not found: ", in_path)

lines <- readLines(in_path, warn = FALSE)

# helper: parse a line like "A 0.1 0.2 ..." into numeric vector
parse_base_line <- function(line, base, L) {
  line <- trimws(line)
  if (!startsWith(line, base)) return(NULL)
  toks <- unlist(strsplit(line, "[\t[:space:]]+"))
  toks <- toks[toks != ""]
  if (length(toks) < (L + 1)) return(NULL)   # base + L numbers
  vals <- suppressWarnings(as.numeric(toks[2:(L+1)]))
  if (any(is.na(vals))) return(NULL)
  vals
}

# softmax over A,C,G,T at each position (column-wise after we build)
softmax_vec <- function(x) {
  x <- x - max(x)
  ex <- exp(x)
  ex / sum(ex)
}

motifs <- list()
i <- 1
while (i <= length(lines)) {
  line <- lines[i]
  if (grepl("^>matrix[0-9]+", line)) {
    # extract name matrix#
    mname <- sub("^>(matrix[0-9]+).*", "\\1", line)

    # expect next 4 non-empty lines are A/C/G/T (sometimes blank lines could exist, so scan forward)
    j <- i + 1
    next_nonempty <- function(k) {
      while (k <= length(lines) && trimws(lines[k]) == "") k <- k + 1
      k
    }
    j <- next_nonempty(j); A <- parse_base_line(lines[j], "A", filter_len)
    j <- next_nonempty(j+1); C <- parse_base_line(lines[j], "C", filter_len)
    j <- next_nonempty(j+1); G <- parse_base_line(lines[j], "G", filter_len)
    j <- next_nonempty(j+1); T <- parse_base_line(lines[j], "T", filter_len)

    if (is.null(A) || is.null(C) || is.null(G) || is.null(T)) {
      stop("Failed to parse A/C/G/T lines for motif ", mname,
           " near line ", i, ". Check file formatting around that motif.")
    }

    # Build a 4 x L matrix in A,C,G,T order, then convert to L x 4 for MEME
    W <- rbind(A, C, G, T)  # 4 x L (log weights)
    # Convert to probabilities per position using softmax across the 4 bases
    P <- apply(W, 2, softmax_vec)  # 4 x L
    P <- t(P)  # L x 4

    motifs[[mname]] <- P
    i <- j + 1
  } else {
    i <- i + 1
  }
}

if (length(motifs) == 0) stop("No motifs parsed. Did not find lines starting with '>matrix'.")

# Write MEME
con <- file(out_path, open = "wt")
on.exit(close(con), add = TRUE)

writeLines("MEME version 4\n", con)
writeLines("ALPHABET= ACGT\n", con)
writeLines("strands: + -\n", con)
writeLines("Background letter frequencies:", con)
writeLines("A 0.25 C 0.25 G 0.25 T 0.25\n", con)

for (nm in names(motifs)) {
  P <- motifs[[nm]]
  writeLines(paste0("MOTIF ", nm), con)
  writeLines(sprintf("letter-probability matrix: alength= 4 w= %d nsites= 20 E= 0", nrow(P)), con)
  for (r in seq_len(nrow(P))) {
    writeLines(paste(sprintf("%.6f", P[r, ]), collapse = " "), con)
  }
  writeLines("", con)
}

cat("DONE ✅\n")
cat("Input:  ", in_path, "\n", sep="")
cat("Output: ", out_path, "\n", sep="")
cat("Motifs parsed: ", length(motifs), "\n", sep="")
cat("Motif width: ", filter_len, "\n", sep="")
cat("Example motifs: ", paste(head(names(motifs), 5), collapse = ", "), "\n", sep="")