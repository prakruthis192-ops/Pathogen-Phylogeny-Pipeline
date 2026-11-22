import os
import subprocess
import matplotlib.pyplot as plt
from Bio import Phylo

CODON_ALIGNED_DNA = "data/processed/h3n2_ha_codon_aligned.fasta"
TREE_PREFIX = "results/trees/h3n2_ha_ml"
TREEFILE = f"{TREE_PREFIX}.treefile"
FIG_OUT = "results/figures/h3n2_phylogeny.png"


def run_iqtree():
    """Infer maximum-likelihood tree with ModelFinder and 1000 ultrafast bootstraps."""
    os.makedirs("results/trees", exist_ok=True)
    print("[IQ-TREE] Running ModelFinder + ML tree inference (1000 ultrafast bootstraps)...")
    cmd = [
        "iqtree",
        "-s", CODON_ALIGNED_DNA,
        "-m", "MFP",
        "-bb", "1000",
        "-nt", "AUTO",
        "-pre", TREE_PREFIX,
        "-redo",
        "-quiet"
    ]
    subprocess.run(cmd, check=True)
    print(f"[IQ-TREE] Tree saved to {TREEFILE}")


def plot_phylogeny():
    """Render a clean publication-style phylogenetic tree with bootstrap support."""
    os.makedirs("results/figures", exist_ok=True)
    tree = Phylo.read(TREEFILE, "newick")
    tree.ladderize(reverse=True)

    fig, ax = plt.subplots(figsize=(12, 8))
    plt.rcParams["font.sans-serif"] = "Arial"

    def get_label(clade):
        if clade.is_terminal():
            return f"  {clade.name}"
        if clade.confidence is not None and clade.confidence >= 70:
            return f"{int(clade.confidence)}"
        return ""

    Phylo.draw(
        tree,
        axes=ax,
        do_show=False,
        label_func=get_label,
        branch_labels=None
    )

    ax.set_title(
        "Maximum-Likelihood Phylogeny of Influenza A (H3N2) Hemagglutinin (HA)\n"
        "ModelFinder Best-Fit + 1,000 Ultrafast Bootstraps",
        fontsize=13,
        fontweight="bold",
        pad=15
    )
    ax.set_xlabel("Substitutions per Site", fontsize=11)
    ax.ticklabel_format(axis="x", style="sci", scilimits=(-2, 2))
    ax.set_ylabel("Isolate", fontsize=11)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", linestyle="--", alpha=0.4)

    for text in ax.texts:
        label_txt = text.get_text().strip()
        if label_txt.isdigit():
            text.set_color("#1f78b4")
            text.set_fontsize(8)
            text.set_fontweight("bold")
        else:
            text.set_fontsize(9)

    plt.tight_layout()
    plt.savefig(FIG_OUT, dpi=300)
    print(f"[Plot] High-resolution phylogeny saved to {FIG_OUT}")


if __name__ == "__main__":
    run_iqtree()
    plot_phylogeny()

