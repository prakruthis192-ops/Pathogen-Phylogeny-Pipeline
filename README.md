# Pathogen Phylogeny & Selection Pipeline

An automated, reproducible bioinformatics workflow for retrieving viral coding sequences (CDS), performing codon-aware multiple sequence alignment, and inferring maximum-likelihood phylogenies.

## Current Progress
- [x] **Day 1:** Reproducible `micromamba` environment (`envs/environment.yml`) and automated NCBI Entrez sequence retrieval (`scripts/fetch_sequences.py`) for Influenza A (H3N2) Hemagglutinin (HA) genes.
- [x] **Day 2:** ORF quality filtering and codon-aware alignment (`scripts/codon_align.py` using `Biopython` + `MAFFT`, yielding a 1,698 bp alignment across 25 intact sequences).
- [ ] **Day 3:** Maximum-likelihood phylogenetic tree inference (`IQ-TREE`) and automated tree visualization.
- [ ] **Day 4:** Evolutionary selection metrics and nucleotide diversity analysis.
- [ ] **Day 5:** Full biological interpretation of clades and end-to-end Bash runner.

## Setup & Usage
```bash
# 1. Create and activate environment
micromamba create -f envs/environment.yml
micromamba activate phylo-pipe

# 2. Fetch raw CDS data & run codon-aware alignment
python3 scripts/fetch_sequences.py
python3 scripts/codon_align.py
```
