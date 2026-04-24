#!/usr/bin/env python3
"""
03_group_by_og.py
-----------------
Reorganizes per-species promoter fasta files into per-OG fasta files
required by promoter_reverse_homology.py.

Reads:
  - OG_MAP_FILTERED (3-col tsv: OG, species, gene_id)
  - Per-species fasta files in BY_SPECIES_DIR

Writes:
  - One fasta file per OG in OG_OUT_DIR
    e.g. OG_OUT_DIR/OG0000001.fa

Sequence matching:
  The gene_id in the map must appear somewhere in the fasta header.
  Headers written by extract_promoters.py look like:
    >Species|gene_id|contig:start-end(strand)|upN|OG=OGxxxxxxx
  We match on the gene_id field (2nd pipe-delimited token).

Usage:
  python 03_group_by_og.py \\
      --map    og_map_filtered.tsv \\
      --in_dir by_species_up800 \\
      --out_dir og_promoters/up800 \\
      --min_seqs 8

  Run twice: once for up800 (used in training) and once for up500
  (used in compute_promoter_rep.py).
"""

import argparse
import os
import sys
from collections import defaultdict


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--map",      required=True,
                    help="3-col filtered OG map (OG, species, gene_id)")
    ap.add_argument("--in_dir",   required=True,
                    help="Directory of per-species fasta files")
    ap.add_argument("--out_dir",  required=True,
                    help="Output directory for per-OG fasta files")
    ap.add_argument("--min_seqs", type=int, default=1,
                    help="Skip OG files with fewer than this many sequences (default: 1)")
    ap.add_argument("--ext",      default=".fa",
                    help="Extension for output files (default: .fa)")
    return ap.parse_args()


def load_og_map(map_path):
    """
    Returns:
      og_to_genes : dict  OG -> list of (species, gene_id)
      gene_to_og  : dict  gene_id -> OG  (for fast lookup)
    """
    og_to_genes = defaultdict(list)
    gene_to_og  = {}
    with open(map_path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 3:
                print(f"[WARN] skipping malformed line {lineno}: {line!r}",
                      file=sys.stderr)
                continue
            og, species, gene_id = parts
            og_to_genes[og].append((species, gene_id))
            gene_to_og[gene_id] = og
    return og_to_genes, gene_to_og


def index_species_fasta(fasta_path):
    """
    Reads a per-species fasta and returns:
      { gene_id : (header, sequence) }
    gene_id is extracted as the 2nd pipe-delimited token in the header,
    which matches the protein_id written by extract_promoters.py.

    Example header:
      >AfumigatusAf293|XP_748023.1|scaffold:1-800(+)|up800|OG=OG0000001
    gene_id would be: XP_748023.1
    """
    index = {}
    header = None
    seq_chunks = []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    gene_id = _extract_gene_id(header)
                    if gene_id:
                        index[gene_id] = (header, "".join(seq_chunks))
                header = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if header is not None:
            gene_id = _extract_gene_id(header)
            if gene_id:
                index[gene_id] = (header, "".join(seq_chunks))
    return index


def _extract_gene_id(header):
    """Extract the 2nd pipe-delimited token from a header string."""
    parts = header.split("|")
    if len(parts) >= 2:
        return parts[1]
    # Fallback: use entire header stripped of whitespace
    return header.split()[0]


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    print(f"Loading OG map: {args.map}", flush=True)
    og_to_genes, gene_to_og = load_og_map(args.map)
    print(f"  {len(og_to_genes):,} OGs, {len(gene_to_og):,} gene entries", flush=True)

    # ── Index all per-species fasta files ─────────────────────
    # Build a combined gene_id -> (header, seq) across all species
    print(f"Indexing per-species fasta files in: {args.in_dir}", flush=True)
    gene_index = {}  # gene_id -> (header, seq)
    fa_files = [f for f in os.listdir(args.in_dir)
                if f.endswith(".fa") or f.endswith(".fasta") or f.endswith(".fa.gz")]
    if not fa_files:
        sys.exit(f"[ERROR] No fasta files found in {args.in_dir}")

    for i, fname in enumerate(sorted(fa_files)):
        fpath = os.path.join(args.in_dir, fname)
        sp_index = index_species_fasta(fpath)
        gene_index.update(sp_index)
        if (i + 1) % 50 == 0:
            print(f"  indexed {i+1}/{len(fa_files)} files, "
                  f"{len(gene_index):,} sequences so far", flush=True)

    print(f"  Total sequences indexed: {len(gene_index):,}", flush=True)

    # ── Write per-OG fasta files ───────────────────────────────
    print(f"Writing per-OG fasta files to: {args.out_dir}", flush=True)
    written   = 0
    skipped   = 0
    not_found = 0

    for og, members in og_to_genes.items():
        out_path = os.path.join(args.out_dir, og + args.ext)
        seqs_written = 0
        with open(out_path, "w") as out:
            for species, gene_id in members:
                if gene_id not in gene_index:
                    not_found += 1
                    continue
                header, seq = gene_index[gene_id]
                if not seq:
                    continue
                out.write(f">{header}\n")
                for j in range(0, len(seq), 60):
                    out.write(seq[j:j+60] + "\n")
                seqs_written += 1

        if seqs_written < args.min_seqs:
            os.remove(out_path)
            skipped += 1
        else:
            written += 1

    print(f"\nSummary:")
    print(f"  OG files written:          {written:,}")
    print(f"  OG files skipped (<{args.min_seqs} seqs): {skipped:,}")
    print(f"  Gene IDs not found:        {not_found:,}")
    if not_found > 0:
        print(f"\n  [WARN] {not_found:,} gene IDs in the map were not found in any fasta.")
        print(  "         Check that --id_attr in extract_promoters.py matches")
        print(  "         the IDs stored in Orthogroups.tsv.")


if __name__ == "__main__":
    main()
