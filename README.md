# Aspergillus
This repository contains trial codes and interim figures for figuring out cis-regulatory networks of Aspergillus.

# 1) Configure paths/parameters
nano ~/bigdata/Aspergillus/scripts/00_config.sh

# 2) Prepare inputs (species lists, 3-col map, filtered map)
bash ~/bigdata/Aspergillus/scripts/01_setup.sh


Scripts:
# Step 1:
First make sure you have these files:
DNA fasta
gff
Orthogroups.tsv

# Step 2:
After that run 02_a_extract_promoters.sh directly to get 800bp promoter sequences for training. Took about 10 mins for Aspergillus and got ~5000000 promoter sequences 

# Step 3:
Run 02_b_extract_promoters_500.sh to get 500bp promoter sequences for final results. Took about 5 mins

# Step 4:
Then follow the steps highlighted in Alan Moses github for cis-regulatory networks. Start with creating conda env for biopython.
I used this script:
conda env list #because i keep forgetting the conda names 
conda create -n promoter_tools python=3.11 biopython -y  #takes a while
conda activate promoter_tools 

conda activate alan_tf_gpu #first install conda based on Alan's instructions

Use this format to run the reverse homology script:
Step 1. Learn a motif-based encoder using reverse homology

$ python promoter_reverse_homology.py test_data/

Step 2. compute the RHIEPA representation for the sequences

$ python compute_promoter_rep.py promoter_rh_interpretableN256L15_f4_b256_final_weights.h5 test_data/

# replace test_data with your input file, in my case I used promoters_500bp_15sp for 15 species

# Step 5:
Run meme.sh job script followed by tomtom_sh and parse_tomtom.sh

# Step 6:
You can use any visualization tools based on your research question. The simplified one I am using is 9_interim_fig.sh

