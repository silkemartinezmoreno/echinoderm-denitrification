#!/usr/bin/env bash
# =============================================================================
# 02. Functional inference with PICRUSt2
# -----------------------------------------------------------------------------
# Manuscript section: Methods, "Functional inference using PICRUSt2".
#
# PICRUSt2 runs from the command line (Linux/macOS, or WSL on Windows). 
# Run this script from the project folder:
#   bash scripts/02_run_picrust2.sh
#
# Inputs (written by script 01):
#   results/sequence_processing/otu_table_picrust2.tsv
#   results/sequence_processing/otu_sequences.fasta
#
# Main output used downstream:
#   results/picrust2/KO_metagenome_out/pred_metagenome_contrib.tsv
# This file is very large (tens of GB) and is not included in the repository.
# Script 03 extracts from it only the denitrification KOs, which are provided
# in data/picrust2/denitrification_KOs_contrib.xlsx.
#
# Installation (once; biom-format is installed together with PICRUSt2):
#   conda create -n picrust2 -c bioconda -c conda-forge picrust2
# =============================================================================

set -euo pipefail

# On Windows, open the Linux environment first with `wsl` and move to the
# project folder, e.g. cd /mnt/c/Users/<user>/Desktop/echinoderm-denitrification
eval "$(conda shell.bash hook)"
conda activate picrust2

seq_dir=results/sequence_processing

# 1. Convert the OTU table to HDF5 BIOM format. For the manuscript this was
# done on Windows in a separate conda environment with biom-format 2.1.17
# (conda create -n biom -c conda-forge -c bioconda biom-format); the same
# command works there or here, inside the PICRUSt2 environment.
biom convert \
  -i "$seq_dir/otu_table_picrust2.tsv" \
  -o "$seq_dir/otu_table_hdf5.biom" \
  --table-type="OTU table" \
  --to-hdf5

# 2. Run PICRUSt2
picrust2_pipeline.py \
  -s "$seq_dir/otu_sequences.fasta" \
  -i "$seq_dir/otu_table_hdf5.biom" \
  -o results/picrust2 \
  -p 20 \
  --stratified \
  --per_sequence_contrib \
  --verbose 2>&1 | tee results/picrust2_run.log

# Script 03 reads the uncompressed contribution file
if [ -f results/picrust2/KO_metagenome_out/pred_metagenome_contrib.tsv.gz ]; then
  gunzip -k results/picrust2/KO_metagenome_out/pred_metagenome_contrib.tsv.gz
fi

echo "PICRUSt2 finished"

# Tip: this run takes many hours. To keep it running after closing the
# terminal, start it inside a screen session:
#   screen -S picrust2
#   bash scripts/02_run_picrust2.sh
