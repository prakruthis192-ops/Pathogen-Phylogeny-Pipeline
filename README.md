# Pathogen Phylogeny & Selection Pipeline

An automated, reproducible bioinformatics workflow for retrieving viral coding sequences (CDS), performing codon-aware multiple sequence alignment, and inferring maximum-likelihood phylogenies.

## Current Progress
- [x] **Day 1:** Reproducible `micromamba` environment (`envs/environment.yml`) and automated NCBI Entrez sequence retrieval (`scripts/fetch_sequences.py`) for Influenza A (H3N2) Hemagglutinin (HA) genes.
- [ ] **Day 2:** Sequence quality filtering and codon-aware alignment (`Biopython` + `MAFFT`).
- [ ] **Day 3:** Maximum-likelihood phylogenetic tree inference (`IQ-TREE`) and automated tree visualization.
- [ ] **Day 4:** Evolutionary selection metrics and nucleotide diversity analysis.
- [ ] **Day 5:** Full biological interpretation of clades and end-to-end Bash runner.

## Setup

