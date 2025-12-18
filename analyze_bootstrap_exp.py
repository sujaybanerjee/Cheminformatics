import numpy as np
import matplotlib.pyplot as plt

# Files generated from run_bootstrap_experiment.py
FILES = {
    "n100": "pvals_n100.txt",
    "n500": "pvals_n500.txt",
    "n1000": "pvals_n1000.txt"
}

def load_pvals(fname):
    """Load p-values from a file into a numpy array."""
    with open(fname) as f:
        vals = [float(line.strip()) for line in f if line.strip() != ""]
    return np.array(vals)


def summarize(vals, label):
    """Print summary statistics for a set of p-values."""
    print(f"\n=== Results for {label} ===")
    print(f"Count: {len(vals)}")
    print(f"Mean: {np.mean(vals):.4f}")
    print(f"Std Dev: {np.std(vals, ddof=1):.4f}")


def plot_hist(vals, label):
    """Save histogram of p-values."""
    plt.figure(figsize=(6, 4))
    plt.hist(vals, bins=20, edgecolor="black")
    plt.title(f"P-value Histogram ({label})")
    plt.xlabel("p-value")
    plt.ylabel("Frequency")
    plt.tight_layout()
    plt.savefig(f"hist_{label}.png", dpi=150)
    plt.close()


if __name__ == "__main__":
    for label, fname in FILES.items():
        vals = load_pvals(fname)
        summarize(vals, label)
        plot_hist(vals, label)

    print("\nSaved histograms:")
    for label in FILES.keys():
        print(f" - hist_{label}.png")