###copyright Alan Moses 2021-2025
##permission to use and copy is granted for research purposes only
##software provided as is with no warrantee of any kind

"""
compute_promoter_rep.py
-----------------------
Encodes per-OG promoter fasta files into phylogenetically-averaged
256-dimensional RHIEPA representations.

Usage:
  python compute_promoter_rep.py weights.h5 og_up500_dir/ \\
      --seq_l 500 --min_seqs 5 --filter_num 256 --filter_len 15 \\
      --out results/rep.csv
"""

import sys
import argparse
import numpy as np
import mol_evol as me
import encoders as my_encoders


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("weights",     help="Path to trained encoder weights (.h5)")
    ap.add_argument("seq_dir",     help="Directory of per-OG fasta files (up500)")
    ap.add_argument("--seq_l",     type=int,   default=500,
                    help="Sequence length for encoding. Paper uses 500 (default: 500). "
                         "Does NOT need to match the training SEQ_L of 800.")
    ap.add_argument("--min_seqs",  type=int,   default=5,
                    help="Minimum orthologs required to compute a representation "
                         "(default: 5). OGs below this are skipped. "
                         "Raise this if you want representations only for "
                         "well-covered gene families.")
    ap.add_argument("--filter_num", type=int,  default=256)
    ap.add_argument("--filter_len", type=int,  default=15)
    ap.add_argument("--out",        default="rep.csv",
                    help="Output CSV path (default: rep.csv)")
    ap.add_argument("--batch_size", type=int,  default=64,
                    help="Batch size for enc.predict_on_batch (default: 64)")
    return ap.parse_args()


def fasta2one_hot(ALN, standardize, seq_l):
    """Convert all sequences in an Alignment to a stacked one-hot array."""
    seqs = []
    for s in ALN.names:
        oh = my_encoders.dna_one_hot2(ALN.seq[s], standardize=seq_l)
        seqs.append(oh)
    return np.array(seqs)


def main():
    args = parse_args()

    # ── Build encoder ──────────────────────────────────────────
    enc = my_encoders.DNAencoder(
        filterN=args.filter_num,
        filterL=args.filter_len,
        trainable=True,
        print_summary=False
    ).interpretable_encoder(seq_len=args.seq_l)

    print(f"Loading weights from: {args.weights}")
    enc.load_weights(args.weights)

    # ── Load OG fasta files ────────────────────────────────────
    print(f"Reading per-OG fasta files from: {args.seq_dir}")
    DB = me.AlignmentDB()
    DB.read_from_fa_dir(args.seq_dir, ext="fa")
    print(f"  {len(DB.names)} OG files loaded")

    # ── Encode each OG by phylogenetic averaging ───────────────
    REP      = {}
    skipped  = 0
    encoded  = 0

    for a in DB.names:
        ALN = DB.data[a]

        # Strip gaps and truncate to last seq_l bp
        # (last seq_l bp = closest to start codon, which is at the 3' end
        #  of the upstream region as extracted by extract_promoters.py)
        clean_names = []
        for s in ALN.names:
            seq = ALN.seq[s].replace("-", "")
            if len(seq) > args.seq_l:
                seq = seq[-args.seq_l:]
            ALN.seq[s] = seq
            if len(seq) >= 10:          # discard near-empty sequences
                clean_names.append(s)
        ALN.names = clean_names

        # Apply minimum orthologs threshold
        # FIX: original code only required >1; we now require min_seqs
        # to avoid noisy representations from poorly-covered OGs.
        if len(ALN.names) < args.min_seqs:
            skipped += 1
            continue

        one_hot = fasta2one_hot(ALN, standardize=args.seq_l, seq_l=args.seq_l)

        # Phylogenetic average: mean over all orthologs
        reps = enc.predict_on_batch(one_hot)
        REP[a] = np.mean(reps, axis=0)
        encoded += 1

    print(f"  Encoded: {encoded}  Skipped (<{args.min_seqs} seqs): {skipped}")

    if not REP:
        sys.exit("[ERROR] No OGs encoded. Check min_seqs threshold and input directory.")

    # ── Write output CSV ───────────────────────────────────────
    dim = len(next(iter(REP.values())))
    print(f"  Representation dimensionality: {dim}")
    print(f"  Writing to: {args.out}")

    with open(args.out, 'w') as out:
        out.write("FILE," + ",".join(str(i) for i in range(dim)) + "\n")
        for a, vec in REP.items():
            out.write(a + "," + ",".join(str(x) for x in vec) + "\n")

    print("Done.")


if __name__ == "__main__":
    main()
