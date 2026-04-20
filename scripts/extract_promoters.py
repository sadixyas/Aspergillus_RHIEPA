#!/usr/bin/env python3
import argparse
import sys
from collections import defaultdict
import re

def parse_attrs(attr_str: str) -> dict:
    d = {}
    for part in attr_str.strip().split(";"):
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            d[k] = v
    return d

def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]

def load_fasta_as_dict(path: str) -> dict:
    seqs = {}
    name = None
    chunks = []
    opener = open
    if path.endswith(".gz"):
        import gzip
        opener = gzip.open
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if name is not None:
            seqs[name] = "".join(chunks)
    return seqs

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", required=True)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--upstream", type=int, default=500)
    ap.add_argument("--species", required=True)
    ap.add_argument("--id_attr", default="protein_id",
                    help="GFF attribute that matches OrthoFinder IDs (fallback is CDS Parent).")
    ap.add_argument("--drop_ambig", action="store_true")
    ap.add_argument("--min_len", type=int, default=50)
    args = ap.parse_args()

    genome = load_fasta_as_dict(args.genome)
    contig_len = {k: len(v) for k, v in genome.items()}

    cds_by_parent = defaultdict(list)         # Parent -> list of CDS
    cds_intervals_by_contig = defaultdict(list)  # contig -> list of (start,end)

    gff_opener = open
    if args.gff.endswith(".gz"):
        import gzip
        gff_opener = gzip.open

    with gff_opener(args.gff, "rt") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            contig, _, ftype, start, end, _, strand, _, attrs = parts
            if ftype != "CDS":
                continue
            if contig not in genome:
                continue
            start = int(start); end = int(end)
            ad = parse_attrs(attrs)
            parent = ad.get("Parent") or ad.get("ID")
            if not parent:
                continue
            for p in parent.split(","):
                cds_by_parent[p].append((contig, start, end, strand, ad))
            cds_intervals_by_contig[contig].append((start, end))

    for contig in cds_intervals_by_contig:
        cds_intervals_by_contig[contig].sort(key=lambda x: (x[0], x[1]))

    def nearest_upstream_cds_end(contig: str, pos: int) -> int:
        best = 0
        for s, e in cds_intervals_by_contig[contig]:
            if e < pos and e > best:
                best = e
        return best

    def nearest_downstream_cds_start(contig: str, pos: int):
        best = None
        for s, _e in cds_intervals_by_contig[contig]:
            if s > pos:
                best = s if best is None else min(best, s)
        return best

    written = 0
    skipped = 0

    with open(args.out, "w") as out:
        for parent, cds_list in cds_by_parent.items():
            contig = cds_list[0][0]
            strand = cds_list[0][3]
            if contig not in contig_len:
                skipped += 1
                continue

            # Start codon proxy from CDS coords
            if strand == "+":
                start_codon = min(x[1] for x in cds_list)
                prom_start = max(1, start_codon - args.upstream)
                prom_end = start_codon - 1
                up_end = nearest_upstream_cds_end(contig, start_codon)
                if up_end > 0:
                    prom_start = max(prom_start, up_end + 1)
            elif strand == "-":
                start_codon = max(x[2] for x in cds_list)
                prom_start = start_codon + 1
                prom_end = min(contig_len[contig], start_codon + args.upstream)
                down_start = nearest_downstream_cds_start(contig, start_codon)
                if down_start is not None:
                    prom_end = min(prom_end, down_start - 1)
            else:
                skipped += 1
                continue

            if prom_end < prom_start:
                skipped += 1
                continue
            if (prom_end - prom_start + 1) < args.min_len:
                skipped += 1
                continue

            # Output ID for matching Orthogroups.tsv
            out_id = parent
            for *_rest, ad in cds_list:
                v = ad.get(args.id_attr)
                if v:
                    out_id = v
                    break

            seq = genome[contig][prom_start-1:prom_end]
            if strand == "-":
                seq = revcomp(seq)
            seq = seq.upper()

            if args.drop_ambig and re.search(r"[^ACGT]", seq):
                skipped += 1
                continue

            header = f">{args.species}|{out_id}|{contig}:{prom_start}-{prom_end}({strand})|up{args.upstream}"
            out.write(header + "\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + "\n")
            written += 1

    print(f"[done] {args.species} wrote={written} skipped={skipped}", file=sys.stderr)

if __name__ == "__main__":
    main()