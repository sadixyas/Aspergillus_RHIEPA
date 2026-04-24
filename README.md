# Inferring Cis-Regulatory Networks in *Aspergillus* from Genome Sequences

This repository implements the RHIEPA (Reverse Homology Interpretable Encoder Phylogenetic Average) pipeline from [Moses et al. 2026 (*Genetics*)](https://doi.org/10.1093/genetics/iyaf209) to infer transcription factor regulatory networks in *Aspergillus* using genome sequences alone — no experimental ChIP-seq or expression data required.

The core idea: orthologous promoter regions across many species share conserved transcription factor binding motifs even when their DNA sequences have diverged too much to align. By training a neural network to recognize orthologs using only their promoter sequences, we can learn biologically meaningful motif-based representations and use them to predict regulatory connections.

---

## Table of Contents

1. [Background](#background)
2. [Requirements](#requirements)
3. [Repository structure](#repository-structure)
4. [Pipeline overview](#pipeline-overview)
5. [Step-by-step instructions](#step-by-step-instructions)
6. [Key parameters and how to set them](#key-parameters-and-how-to-set-them)
7. [Expected outputs and how to interpret them](#expected-outputs-and-how-to-interpret-them)
8. [Troubleshooting](#troubleshooting)
9. [Citation](#citation)

---

## Background

### What is a promoter region?
In fungi, the promoter region is the DNA sequence immediately upstream of a gene's start codon. It contains binding sites for transcription factors (TFs) — proteins that control when and how much a gene is expressed. This pipeline extracts up to 800 bp upstream of each gene as a proxy for its promoter.

### What is a transcription factor motif (PWM)?
A position weight matrix (PWM) describes the sequence preference of a TF — which nucleotides it prefers at each position of its binding site. This pipeline learns 256 PWMs directly from sequence data without needing prior knowledge of which TFs exist in *Aspergillus*.

### What is reverse homology / contrastive learning?
Contrastive learning trains a model to make similar representations for related sequences (orthologs) and different representations for unrelated sequences. Here, "related" means "orthologous promoters across species": they regulate the same gene and are therefore expected to share functional regulatory elements. This signal is used to learn which sequence features (motifs) are conserved and therefore likely functional.

### What is phylogenetic averaging (RHIEPA)?
Once the encoder is trained, we encode the promoter from every orthologous species and take the average across all of them. This averaging amplifies conserved motif signals and suppresses noise from non-functional motif matches that occur by chance in any single genome, addressing the classic "futility theorem" problem of motif scanning in single genomes.

---

### Software

**Conda environments** — set up before running anything:

```bash
# Environment 1: promoter extraction tools
conda create -n promoter_tools python=3.11 biopython -y
conda activate promoter_tools

# Environment 2: TensorFlow + GPU for model training
# Follow Alan Moses lab instructions at:
# https://github.com/alanmoses/promoter_reverse_homology
conda activate alan_tf_gpu

# Environment 3: MEME Suite for TomTom motif comparison
conda create -n meme_suite -c bioconda meme -y
```

**R packages** — run inside R before step 6:
```r
install.packages(c("readr", "dplyr", "stringr"))
```

**Required input files per species:**
- Genome DNA fasta (`.fa` or `.fa.gz`)
- Genome annotation in GFF3 format (`.gff3` or `.gff`)
- `Orthogroups.tsv` from OrthoFinder, one row per orthogroup, one column per species

**Motif databases** — download before steps 7 and 8:
- JASPAR 2026 fungi (MEME format): https://jaspar.elixir.no/downloads/
- Cis-BP 2.0 for *Aspergillus fumigatus*: https://cisbp.ccbr.utoronto.ca
  - Select species, download "PWMs (MEME)" format
  - Save to `motif_db/cisbp2_aspergillus.meme`

---

## Repository structure

```
scripts/
├── 00_config.sh                  <- All paths and parameters live here
├── 01_setup.sh                   <- Filter orthogroups, build OG map
├── 02_extract_promoters.sh       <- Extract per-species promoter fastas
├── 03_group_by_og.py             <- Reorganise into per-OG fastas
├── 03_group_by_og.sh             <- SLURM wrapper for above
├── 04_train_rh.sh                <- Train the interpretable encoder (GPU)
├── 05_compute_rep.sh             <- Compute RHIEPA representations
├── 06_meme.sh                    <- Convert learned motifs to MEME format
├── 07_tomtom.sh                  <- Compare motifs to JASPAR + Cis-BP
├── 08_parse_tomtom.sh            <- Parse results with correct thresholds
├── extract_promoters.py          <- Promoter extraction from GFF + genome
├── promoter_reverse_homology.py  <- Contrastive learning training script
├── encoders.py                   <- Neural network architecture (PWM encoder)
├── compute_promoter_rep.py       <- RHIEPA representation computation
├── mol_evol.py                   <- Sequence I/O utilities (Moses lab)
├── matrixlayer_to_meme.R         <- Convert learned PWMs to MEME format
└── parse_tomtom.R                <- Parse TomTom output with correct cutoffs

motif_db/
├── JASPAR2026_CORE_fungi_non-redundant_pfms_meme.txt
└── cisbp2_aspergillus.meme       <- Download separately
```

---

## Pipeline overview

```
Genome fastas + GFF files
        |
        v
[01] Filter orthogroups (01_setup.sh)
        |  Orthogroups.tsv -> filtered OG map
        v
[02] Extract promoter sequences (02_extract_promoters.sh)
        |  Per-species: by_species_up800/, by_species_up500/
        v
[03] Group sequences by orthogroup (03_group_by_og.sh)
        |  Per-OG: og_promoters/up800/, og_promoters/up500/
        v
[04] Train interpretable encoder (04_train_rh.sh)  <- GPU, 1-7 days
        |  Learns 256 PWMs via contrastive learning on up800 sequences
        v
[05] Compute RHIEPA representations (05_compute_rep.sh)
        |  Encodes up500 sequences -> phylogenetically averaged 256-dim vectors
        v
[06] Convert learned motifs to MEME format (06_meme.sh)
        |  matrix_layer.txt -> motifs.meme
        v
[07] Compare motifs to known TF databases (07_tomtom.sh)
        |  TomTom search vs JASPAR fungi + Cis-BP 2.0
        v
[08] Parse and threshold results (08_parse_tomtom.sh)
        |  Strict: E-value < 4e-3 (for counting matched motifs)
        |  Relaxed: q-value <= 0.05 (for heatmap TF name labels)
        v
Downstream analysis (UMAP, heatmaps, GO enrichment, etc.)
```

---

## Step-by-step instructions

### Step 0: Configure all paths and parameters

**Edit `00_config.sh` first.** This is the single source of truth for the entire pipeline — every other script sources it. You should only need to edit this one file to adapt the pipeline to your system.

```bash
nano scripts/00_config.sh
```

Key things to set:
- `ROOT`: your project root directory
- `THRESH`: minimum number of species an orthogroup must be present in.
  With 406 Aspergillus species, I used **264** (~65% coverage).
  Do not set this too low, poorly covered OGs provide weak training signal.
- `DB_CISBP`: path to your downloaded Cis-BP MEME file
- Conda environment names if they differ from the defaults

> **Why 65% coverage?** OGs present in fewer species give the contrastive
> learner fewer examples of what is conserved. Very low-coverage OGs
> (present in fewer than 20% of species) are mostly noise and should be
> excluded.

---

### Step 1: Filter orthogroups

```bash
sbatch scripts/01_setup.sh
# Runtime: ~10 minutes
# Output:  $ORTHO_DIR/og_map_filtered.tsv
```

This script reads `Orthogroups.tsv` from OrthoFinder and:
1. Lists all species in the dataset
2. Selects orthogroups present in at least `THRESH` species
3. Takes the first gene when a species has paralogs (does **not** require single-copy)
4. Writes a 3-column map: `OG | species | gene_id`

The diagnostic output will tell you how many OGs pass the filter. For Aspergillus with `THRESH=264`, I had roughly 3,000–8,000 OGs.

---

### Step 2: Extract promoter sequences

```bash
sbatch scripts/02_extract_promoters.sh
# Runtime: ~10-20 minutes for 400 species
# Output:  $BY_SPECIES_UP800/ and $BY_SPECIES_UP500/
#          One fasta file per species
```

For each species, extracts up to 800 bp (or 500 bp) upstream of every annotated gene start codon, clipped at the nearest upstream coding sequence to avoid overlap with adjacent genes. Sequences containing ambiguous bases (Ns) are discarded.

> **Note on two lengths:** We extract both 800 bp and 500 bp. The 800 bp
> sequences are used during training (more context for learning motifs).
> The 500 bp sequences are used for final representations (reduces
> contamination from adjacent gene UTRs).

---

### Step 3: Group sequences by orthogroup

```bash
sbatch scripts/03_group_by_og.sh
# Runtime: ~30-60 minutes
# Output:  $OG_UP800/ and $OG_UP500/
#          One fasta file per OG, e.g. OG0000001.fa
```

The training script expects one fasta file per orthogroup where all sequences in the file are orthologs of each other. This step reorganises the per-species files from Step 2 into that format.

If you see many "Gene IDs not found" warnings, there is likely a mismatch between the IDs in `Orthogroups.tsv` and the IDs in your GFF files. Check that `--id_attr` in `02_extract_promoters.sh` matches the attribute OrthoFinder used to build the orthogroup table (commonly `protein_id`, `Name`, or `ID`).

---

### Step 4: Train the interpretable encoder

```bash
sbatch scripts/04_train_rh.sh
# Runtime: 1-7 days depending on GPU and number of OGs
# Output:  *_final_weights.h5
#          *_matrix_layer.txt
#          *_training_info.csv
```

This is the core machine learning step. The encoder learns 256 position weight matrices by training to distinguish orthologous promoter groups from unrelated sequences. Each epoch processes all sequences once.

Monitor training by checking the loss in the training CSV:
```bash
grep "^loss" *_training_info.csv
```
Loss should decrease from ~5.5 toward ~2.5–3.5. If it plateaus above 4.5, the model may not have enough data or the OG files may have quality issues.

> **GPU memory:** If you get out-of-memory errors, reduce `batch_size`
> in `promoter_reverse_homology.py` from 256 to 128.

---

### Step 5: Compute RHIEPA representations

```bash
sbatch scripts/05_compute_rep.sh
# Runtime: ~1-3 hours (CPU, no GPU needed)
# Output:  results/rep.csv
#          Rows = OGs, columns = 256 motif signals (phylogenetically averaged)
```

For each OG, encodes all orthologous 500 bp promoter sequences and takes their mean across species. The output `rep.csv` is used for all downstream analysis.

---

### Step 6: Convert learned motifs to MEME format

```bash
sbatch scripts/06_meme.sh
# Runtime: < 5 minutes
# Output:  results/motifs.meme
```

Verify the output has 256 motifs:
```bash
grep -c "^MOTIF" results/motifs.meme
# Should print: 256
```

---

### Step 7: Compare learned motifs to known TF databases

```bash
sbatch scripts/07_tomtom.sh
# Runtime: ~30 minutes
# Output:  results/tomtom_jaspar/tomtom.tsv
#          results/tomtom_cisbp/tomtom.tsv (if Cis-BP file is present)
```

Searches both JASPAR fungi and Cis-BP 2.0. Both databases are needed because different TF binding specificities are deposited in different databases. Searching only JASPAR will systematically undercount matched motifs for Aspergillus.

---

### Step 8: Parse TomTom results

```bash
sbatch scripts/08_parse_tomtom.sh
# Runtime: < 5 minutes
# Output:  results/parsed_tomtom/matched_strict_E4e3.tsv
#          results/parsed_tomtom/matched_relaxed_q05.tsv
#          results/parsed_tomtom/summary.txt
```

Applies significance thresholds:
- **E-value < 4e-3**: for counting how many of the 256 motifs match a known TF
- **q-value <= 0.05**: for assigning TF names to motifs in heatmaps

The `summary.txt` file reports the final count directly.

---

## Key parameters and how to set them

| Parameter | File | Default | Guidance |
|-----------|------|---------|----------|
| `THRESH` | `00_config.sh` | 264 | ~65% of species count. Raise for stricter coverage. |
| `FILTER_NUM` | `00_config.sh` | 256 | Number of PWMs. Match paper; changing requires retraining. |
| `FILTER_LEN` | `00_config.sh` | 15 | PWM width in bp. Match paper. |
| `FAM_SIZE` | `00_config.sh` | 4 | Orthologs averaged per contrastive step. |
| `TRAIN_SEQ_L` | `00_config.sh` | 800 | Sequence length during training. |
| `REP_SEQ_L` | `00_config.sh` | 500 | Sequence length for representations. Keep at 500. |
| `MIN_ORTHOLOGS` | `00_config.sh` | 5 | Min orthologs to compute a representation. |
| `EVALUE_STRICT` | `00_config.sh` | 0.004 | For counting matched motifs. Do not change. |
| `QVALUE_RELAXED` | `00_config.sh` | 0.05 | For heatmap TF labels. Do not change. |

---

## Expected outputs and how to interpret them

### How many motifs should match known TFs?
With ~400 Aspergillus species and both databases searched, expect **30–60** of the 256 learned PWMs to match known motifs at E-value < 4e-3. The paper found 63 for *S. cerevisiae* with 89 species. Fewer species and less complete Cis-BP coverage for Aspergillus means a somewhat lower number is expected and acceptable.

### What does the UMAP show?
After computing `rep.csv`, project the 256-dimensional representations to 2D using UMAP. Promoters of co-regulated genes (ribosome biogenesis, cell cycle, amino acid biosynthesis) should cluster together. A uniform UMAP with no clear clusters indicates the model did not converge — check the training loss curve.

### What does motif collapse look like?
If one or two TF IDs dominate the TomTom bar plot with far more matches than others, the model has collapsed — many of the 256 PWMs converged to the same sequence preference. This is caused by BatchNorm being applied in the wrong position in the encoder. The corrected `encoders.py` in this repository fixes this.

---

## Troubleshooting

**Training loss does not decrease below 4.5**
Check how many OG files have at least 8 sequences:
```bash
for f in og_promoters/up800/*.fa; do grep -c "^>" "$f"; done | awk '$1>=8' | wc -l
```
If fewer than 1,000 OGs pass this threshold, your `THRESH` setting may be too high or the single-copy filter removed too many families.

**TomTom returns very few significant hits**
Confirm you have both JASPAR and Cis-BP databases and that `07_tomtom.sh` ran against both. Check `summary.txt` : if there are many hits at q<=0.05 but few at E<4e-3, motifs are marginal and the model needs retraining with the corrected encoder.

**"Gene IDs not found" warnings in Step 3**
The gene IDs in `Orthogroups.tsv` do not match the `protein_id` attribute in your GFF files. Try `--id_attr Name` or `--id_attr ID` in `02_extract_promoters.sh`.

**CUDA out of memory**
Reduce `batch_size` from 256 to 128 in `promoter_reverse_homology.py`.

**OG fasta files have only 1–2 sequences**
Your `THRESH` is too low — sparse OGs are being included. Raise `THRESH` to ~65% of your species count and rerun from Step 1.

---

## Citation

Original pipeline code:
- https://github.com/alanmoses/promoter_reverse_homology
- https://github.com/stajichlab/Y1000_orthologs
