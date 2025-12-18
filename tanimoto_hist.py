import sys
import pandas as pd
import matplotlib.pyplot as plt


def make_histograms(tanimoto_output_file, sunet_id="sujayb23"):
    # Load the pairwise results
    df = pd.read_csv(
        tanimoto_output_file,
        header=None,
        names=["drugA", "drugB", "tanimoto", "share"]
    )

    # Extract the three distributions
    all_vals = df["tanimoto"].values
    shared_vals = df[df["share"] == 1]["tanimoto"].values
    notshared_vals = df[df["share"] == 0]["tanimoto"].values

    # 1. All tanimoto values
    plt.figure(figsize=(8,6))
    plt.hist(all_vals, bins=40, edgecolor="black")
    plt.title(f"{sunet_id} All")
    plt.xlabel("Tanimoto Similarity")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig("all_tanimoto.png")
    plt.close()

    # 2. Shared target
    plt.figure(figsize=(8,6))
    plt.hist(shared_vals, bins=40, edgecolor="black", color="green")
    plt.title(f"{sunet_id} Shared")
    plt.xlabel("Tanimoto Similarity")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig("shared_tanimoto.png")
    plt.close()

    # 3. Not shared
    plt.figure(figsize=(8,6))
    plt.hist(notshared_vals, bins=40, edgecolor="black", color="red")
    plt.title(f"{sunet_id} Not Shared")
    plt.xlabel("Tanimoto Similarity")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig("notshared_tanimoto.png")
    plt.close()


if __name__ == "__main__":
    tanimoto_output_file = sys.argv[1]
    make_histograms(tanimoto_output_file)

