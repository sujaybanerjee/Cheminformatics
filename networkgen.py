import sys
import argparse
import pandas as pd
from chemoUtils import load_drugs_fpts, load_protein_ligands
from pvalue import bootstrap_pvalue 

def load_network_proteins(protein_nodes_file):
    """
    Load list of uniprot_accession IDs from protein_nodes.csv
    Input: protein_nodes.csv : uniprot_accession, uniprot_id, indications
    Returns:
        list of uniprot_accesions
    """
    uniprots = pd.read_csv(protein_nodes_file)
    return list(uniprots["uniprot_accession"])

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", type=int, default=500, help="num of bootstrap iter")
    parser.add_argument("-r", type=int, default=214, help="seed")
    parser.add_argument("drugs_file")
    parser.add_argument("targets_file")
    parser.add_argument("protein_nodes_file")
    args = parser.parse_args()

    #load data
    drug_to_fpt = load_drugs_fpts(args.drugs_file)
    protein_to_ligands = load_protein_ligands(args.targets_file)
    all_drug_ids = list(drug_to_fpt)
    network_proteins = load_network_proteins(args.protein_nodes_file)

    out_file = "network_edgelist.txt"
    out = open(out_file, "w")
    #p-vals for each unique pair
    for i in range(len(network_proteins)):
        for j in range(i + 1,len(network_proteins)):  
            proteinA = network_proteins[i]
            proteinB = network_proteins[j]
            pval = bootstrap_pvalue(proteinA,proteinB,drug_to_fpt, protein_to_ligands,all_drug_ids,num_iter=args.n, seed=args.r)
            if pval<=0.05:
                out.write(f"{proteinA} {proteinB}\n")
    out.close()