
import sys
import pandas as pd
from chemoUtils import load_drugs_fpts, load_protein_ligands, tanimoto


def compute_tanimoto_output(drugs_file, targets_file, output_file):
    """
    Computes pairwise Tanimoto similarities for all unique drug pairs
    Inputs:
        drugs_file: drugs.csv:
            db_id
            maccs 
        targets_file: targets.csv:
            db_id
            uniprot_accession
        output_file: output CSV file
            drugA,drugB,tanimoto_score,share_target(0/1)
    """
    #load data
    drug_to_fp = load_drugs_fpts(drugs_file)
    protein_to_drugs = load_protein_ligands(targets_file)

    drug_to_proteins = {}
    for protein, drugset in protein_to_drugs.items():
        for drug in drugset:
            if drug not in drug_to_proteins:
                drug_to_proteins[drug] = set()
            drug_to_proteins[drug].add(protein)
    all_drugs = list(drug_to_fp.keys())

    #compute all unique drug pairs
    with open(output_file,"w") as out:
        for i in range(len(all_drugs)):
            for j in range(i+1, len(all_drugs)):
                drugA = all_drugs[i]
                drugB = all_drugs[j]
                #tanimoto
                score = tanimoto(drug_to_fp[drugA], drug_to_fp[drugB])
                #shared target
                protsA = drug_to_proteins.get(drugA, set())
                protsB = drug_to_proteins.get(drugB, set())
                share_target = 1 if (protsA & protsB) else 0

                out.write(f"{drugA},{drugB},{score:.6f},{share_target}\n")

if __name__ == "__main__":
    drugs_file = sys.argv[1]
    targets_file = sys.argv[2]
    output_file = sys.argv[3]

    compute_tanimoto_output(drugs_file, targets_file, output_file)