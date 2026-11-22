import os
from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO

ALIGNED_AA = "data/processed/h3n2_ha_aligned_aa.fasta"
CSV_OUT = "results/mutation_summary.csv"
FIG_OUT = "results/figures/ha_mutation_profile.png"


def find_mutations():
    """Compare each aligned protein sequence against the cohort consensus."""
    records = list(SeqIO.parse(ALIGNED_AA, "fasta"))
    seq_len = len(records[0].seq)

    # 1. Build consensus sequence (most frequent amino acid at each position)
    consensus = []
    for pos in range(seq_len):
        column_residues = [str(r.seq)[pos] for r in records]
        most_common_aa = Counter(column_residues).most_common(1)[0][0]
        consensus.append(most_common_aa)
    consensus_str = "".join(consensus)

    # 2. Detect point mutations relative to consensus (1-indexed protein positions)
    mutation_rows = []
    for rec in records:
        seq_str = str(rec.seq)
        for pos in range(seq_len):
            ref_aa = consensus_str[pos]
            alt_aa = seq_str[pos]
            if alt_aa != ref_aa and alt_aa != "-":
                mutation_rows.append({
                    "Isolate": rec.id,
                    "Position": pos + 1,
                    "Ref_AA": ref_aa,
                    "Alt_AA": alt_aa,
                    "Mutation": f"{ref_aa}{pos + 1}{alt_aa}"
                })

    df = pd.DataFrame(mutation_rows)
    os.makedirs("results", exist_ok=True)
    df.to_csv(CSV_OUT, index=False)
    print(f"[Analysis] Found {len(df)} total non-consensus substitutions across {len(records)} isolates.")
    print(f"[Saved] Detailed mutation table saved to {CSV_OUT}")
    return df, len(records)


def plot_mutations(df, total_isolates):
    """Plot mutation counts per isolate and top recurring substitutions."""
    os.makedirs("results/figures", exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    # Left panel: Top recurring amino acid mutations
    top_muts = df["Mutation"].value_counts()
    ax1.barh(top_muts.index[::-1], top_muts.values[::-1], color="#2b5c8f")
    ax1.set_title("Recurring HA Amino Acid Substitutions", fontweight="bold")
    ax1.set_xlabel(f"Number of Isolates (out of {total_isolates})")
    ax1.grid(axis="x", linestyle="--", alpha=0.4)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    # Right panel: Mutations per isolate
    per_isolate = df["Isolate"].value_counts()
    ax2.barh(per_isolate.index[::-1], per_isolate.values[::-1], color="#d95f02")
    ax2.set_title("Divergent Isolates (Mutations vs. Consensus)", fontweight="bold")
    ax2.set_xlabel("Number of Amino Acid Differences")
    ax2.tick_params(axis="y", labelsize=8)
    ax2.grid(axis="x", linestyle="--", alpha=0.4)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig(FIG_OUT, dpi=300)
    print(f"[Plot] Mutation profile saved to {FIG_OUT}")


if __name__ == "__main__":
    mut_df, n_seqs = find_mutations()
    plot_mutations(mut_df, n_seqs)

