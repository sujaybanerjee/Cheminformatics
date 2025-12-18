import subprocess
import random
import sys

# CONFIG ----------------------------
NUM_SEEDS = 100

PROTEIN_A = "P54577"
PROTEIN_B = "Q7RTX0"

# Use YOUR actual paths here
DRUGS = "data/drugs.csv"
TARGETS = "data/targets.csv"

ITER_COUNTS = [100, 500, 1000]
OUTPUT_FILES = {
    100:  "pvals_n100.txt",
    500:  "pvals_n500.txt",
    1000: "pvals_n1000.txt"
}
# -----------------------------------

# Generate 100 distinct random seeds
seeds = random.sample(range(1_000_000), NUM_SEEDS)

python_cmd = sys.executable  # calls the same python you're using now

print("Running 300 p-value computations...")
print("-----------------------------------")

for N in ITER_COUNTS:
    outfile = OUTPUT_FILES[N]
    with open(outfile, "w") as f:
        print(f"\nRunning n = {N} iterations...")
        for s in seeds:
            cmd = [
                python_cmd, "pvalue.py",
                "-n", str(N),
                "-r", str(s),
                DRUGS, TARGETS,
                PROTEIN_A, PROTEIN_B
            ]

            result = subprocess.run(cmd, capture_output=True, text=True)

            if result.returncode != 0:
                # Something went wrong inside pvalue.py
                print(f"[ERROR] seed={s} n={N}")
                print(result.stderr.strip())
                # Don't write anything to file for this failed run
                continue

            pval_str = result.stdout.strip()

            if not pval_str:
                # pvalue.py ran but printed nothing
                print(f"[WARN] Empty output for seed={s}, n={N}")
                continue

            f.write(pval_str + "\n")
            print(f"seed={s} n={N} → p={pval_str}")

print("\nDone! Files generated:")
for N, fname in OUTPUT_FILES.items():
    print(f"- {fname}")