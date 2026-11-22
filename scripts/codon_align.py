import os
import re
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

RAW_FASTA = "data/raw/h3n2_ha_raw.fasta"
CLEAN_DNA = "data/processed/h3n2_ha_clean_dna.fasta"
CLEAN_AA = "data/processed/h3n2_ha_clean_aa.fasta"
ALIGNED_AA = "data/processed/h3n2_ha_aligned_aa.fasta"
CODON_ALIGNED_DNA = "data/processed/h3n2_ha_codon_aligned.fasta"


def parse_strain_name(record):
    """Extract Accession + Strain/Year from messy NCBI headers."""
    acc = record.id.split(".")[0]
    match = re.search(r"\((A/[^)]+)\)", record.description)
    if match:
        strain = re.sub(r"[^A-Za-z0-9_.-]", "_", match.group(1))
        return f"{acc}_{strain}"
    return acc


def filter_and_translate():
    """Keep intact coding sequences (divisible by 3, no internal stop codons)."""
    os.makedirs("data/processed", exist_ok=True)
    valid_dna = []
    valid_aa = []
    seen_ids = set()

    for rec in SeqIO.parse(RAW_FASTA, "fasta"):
        seq_str = str(rec.seq).upper()

        # Skip ambiguous bases or non-triplet lengths
        if "N" in seq_str or len(seq_str) % 3 != 0:
            continue

        # Trim terminal stop codon if present
        if seq_str[-3:] in {"TAA", "TAG", "TGA"}:
            seq_str = seq_str[:-3]

        dna_seq = Seq(seq_str)
        aa_seq = dna_seq.translate()

        # Skip if internal stop codons exist
        if "*" in str(aa_seq):
            continue

        clean_id = parse_strain_name(rec)
        if clean_id in seen_ids:
            continue
        seen_ids.add(clean_id)

        valid_dna.append(SeqRecord(dna_seq, id=clean_id, description=""))
        valid_aa.append(SeqRecord(aa_seq, id=clean_id, description=""))

    SeqIO.write(valid_dna, CLEAN_DNA, "fasta")
    SeqIO.write(valid_aa, CLEAN_AA, "fasta")
    print(f"[QC] Retained {len(valid_dna)} intact coding sequences.")
    return {r.id: str(r.seq) for r in valid_dna}


def run_mafft_protein():
    """Align translated amino acid sequences using MAFFT."""
    print("[MAFFT] Aligning translated amino acid sequences...")
    cmd = ["mafft", "--auto", "--quiet", CLEAN_AA]
    with open(ALIGNED_AA, "w") as out_f:
        subprocess.run(cmd, stdout=out_f, check=True)


def back_translate_to_codons(dna_dict):
    """Map aligned protein gaps back onto original nucleotide triplets."""
    codon_records = []
    for aa_rec in SeqIO.parse(ALIGNED_AA, "fasta"):
        raw_dna = dna_dict[aa_rec.id]
        codon_aligned = []
        codon_idx = 0

        for aa_char in str(aa_rec.seq):
            if aa_char == "-":
                codon_aligned.append("---")
            else:
                codon_aligned.append(raw_dna[codon_idx : codon_idx + 3])
                codon_idx += 3

        codon_records.append(
            SeqRecord(Seq("".join(codon_aligned)), id=aa_rec.id, description="")
        )

    SeqIO.write(codon_records, CODON_ALIGNED_DNA, "fasta")
    print(f"[Back-Translate] Saved codon alignment ({len(codon_records[0].seq)} bp) to {CODON_ALIGNED_DNA}")


if __name__ == "__main__":
    dna_map = filter_and_translate()
    run_mafft_protein()
    back_translate_to_codons(dna_map)

