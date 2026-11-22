#!/usr/bin/env bash
set -euo pipefail

echo "=================================================="
echo "  Pathogen Phylogeny & Selection Pipeline (H3N2)"
echo "=================================================="

echo "[1/4] Retrieving H3N2 Hemagglutinin CDS records from NCBI..."
python3 scripts/fetch_sequences.py

echo "[2/4] Filtering ORFs & building codon-aware alignment (MAFFT)..."
python3 scripts/codon_align.py

echo "[3/4] Inferring Maximum-Likelihood phylogeny (IQ-TREE)..."
python3 scripts/build_and_plot_tree.py

echo "[4/4] Profiling amino acid substitutions vs. cohort consensus..."
python3 scripts/mutation_analysis.py

echo "=================================================="
echo "  Pipeline finished! Outputs saved in results/"
echo "=================================================="
