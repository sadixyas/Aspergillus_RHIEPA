# Inferring Cis-Regulatory Networks in *Aspergillus*

Adapts the RHIEPA pipeline ([Moses et al. 2026](https://doi.org/10.1093/genetics/iyaf209)) to infer TF regulatory networks in *Aspergillus* from genome sequences alone. The model learns which promoter motifs are conserved across hundreds of species and uses that to predict which transcription factors regulate which genes.

---

## Before you start

**Input files required:**
- Genome fasta + GFF annotation for each species
- `Orthogroups.tsv` from OrthoFinder
- Species panel file (one proteome basename per line, no header, no reference species)
- JASPAR 2026 fungi MEME file → https://jaspar.elixir.no/downloads/
- Cis-BP 2.0 for *A. fumigatus* MEME file → https://cisbp.ccbr.utoronto.ca

**Conda environments:**
```bash
conda create -n promoter_tools python=3.11 biopython -y
conda env create -f alan_tf_gpu.yaml        # for training + GPU
conda create -n meme_suite -c bioconda meme -y
```

**R packages:**
```bash
module load R/4.5.0
Rscript -e 'install.packages(c("readr","dplyr","stringr"))'
```

**All scripts directory must contain these files:**
`mol_evol.py`, `tf_tf.py`, `encoders.py`, `extract_promoters.py`,
`03_group_by_og.py`, `promoter_reverse_homology.py`, `compute_promoter_rep.py`,
`matrixlayer_to_meme.R`, `parse_tomtom.R`

---

## Configuration

**Edit `00_config.sh` before running anything.** Every script sources this file, it is the only place you need to set paths and parameters.

Key settings to update:
| Parameter | What to set |
|-----------|-------------|
| `ROOT` | Your project root directory |
| `THRESH` | Min species per orthogroup. With 406 species it is set to **264** (~65%) |
| `DB_CISBP` | Path to your Cis-BP MEME file |
| `SCRIPTS_DIR` | Directory containing all `.py` and `.R` files |


---

## Pipeline

### Step 1: Filter orthogroups
```bash
sbatch scripts/01_setup.sh
```
Reads `Orthogroups.tsv`, keeps OGs present in at least `THRESH` species, and
writes a 3-column map (`OG | species | gene_id`). Takes ~10 min.

### Step 2: Extract promoter sequences
```bash
sbatch scripts/02_extract_promoters.sh
```
Extracts up to 800 bp and 500 bp upstream of each gene for every species.
Sequences with ambiguous bases (Ns) are discarded. Takes ~20 min.

### Step 3: Group sequences by orthogroup
```bash
sbatch scripts/03_group_by_og.sh
```
Reorganises per-species fastas into per-OG fastas required by the training
script. Produces `og_promoters/up800/` and `og_promoters/up500/`. Takes ~1 hour.

### Step 4: Train the encoder
```bash
sbatch scripts/04_train_rh.sh
```
Trains a neural network with 256 learnable PWMs using contrastive learning on
the 800 bp sequences. Runs up to 100 epochs. Takes 1-7 days depending on GPU
and OG count. Monitor training loss in `*_training_info.csv`, it should
decrease from ~5.5 to ~2.5–3.5.

### Step 5: Compute representations
```bash
sbatch scripts/05_compute_rep.sh
```
Encodes 500 bp promoters using the trained model and averages across all orthologs
(phylogenetic averaging). Writes `results/rep.csv`, one row per OG, 256 columns.
This file is the input to all downstream analysis. Takes ~2 hours (CPU only).

### Step 6: Convert motifs to MEME format
```bash
sbatch scripts/06_meme.sh
```
Converts the 256 learned PWMs to MEME format for TomTom. Verify:
```bash
grep -c "^MOTIF" results/motifs.meme   # should print 256
```

### Step 7: Compare motifs to TF databases
```bash
sbatch scripts/07_tomtom.sh
```
Runs TomTom against both JASPAR fungi and Cis-BP 2.0. Both databases are
required (many *Aspergillus* TF specificities are only in Cis-BP).

### Step 8: Parse results
```bash
sbatch scripts/08_parse_tomtom.sh
```
- **E-value < 4e-3** — for counting matched motifs 
- **q-value ≤ 0.05** — for assigning TF names in heatmaps

Read the summary:
```bash
cat results/parsed_tomtom/summary.txt
```

---

## Troubleshooting

**Training loss stays above 4.5**
Check how many OG files have ≥ 8 sequences — if fewer than 1,000, your OG
files may not have been built correctly or THRESH is too high.
```bash
for f in $OG_UP800/*.fa; do grep -c "^>" "$f"; done | awk '$1>=8' | wc -l
```

**"Gene IDs not found" warnings in step 3**
The IDs in `Orthogroups.tsv` don't match those in your GFF files. Try
changing `--id_attr protein_id` to `--id_attr Name` or `--id_attr ID`
in `02_extract_promoters.sh`.

**GPU out of memory in step 4**
Reduce `batch_size` from 256 to 128 in `promoter_reverse_homology.py`.

**Fewer than 30 motif matches after retraining**
Confirm Cis-BP ran: `wc -l results/tomtom_cisbp/tomtom.tsv`. If empty,
the path in `DB_CISBP` in `00_config.sh` is wrong.

---

## Citation

> Moses AM, Stajich JE, Gasch AP, Knowles DA. (2026). Inferring fungal cis-regulatory networks from genome sequences via unsupervised and interpretable representation learning. *Genetics*, 232(1), iyaf209.
