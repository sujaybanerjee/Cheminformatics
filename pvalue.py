import sys
import argparse
import random
from chemoUtils import load_drugs_fpts, load_protein_ligands, tanimoto


def compute_tsummary(ligand_ids_A, ligand_ids_B, drug_to_fpt, cutoff=0.5):
    """
    Compute Tsummary between two ligands by summing tanimoto above cutoff
    Inputs:
        ligand_ids_A: list of drug IDs that bind protein A
        ligand_ids_B: list of drug IDs that bind protein B
        drug_to_fpt: dict that maps drug ID to its fingerprint set
        cutoff: float val where only tanimoto values greater than this are included
    Returns:
        Tsummary val: sum of Tc > cutoff
    """
    sum = 0.0
    for drugA in ligand_ids_A:
        fptA = drug_to_fpt.get(drugA)
        for drugB in ligand_ids_B:
            fptB = drug_to_fpt.get(drugB)
            Tc = tanimoto(fptA, fptB)
            if Tc > cutoff:
                sum += Tc
    return sum


def bootstrap_pvalue(proteinA_id,proteinB_id,drug_to_fpt,protein_to_ligands, all_drug_ids, num_iter,seed):
    """
    Compute bootstrap p-value for ligand-set similarity between two proteins
    Inputs:
        proteinA_id: uniprot accession ID for protein A
        proteinB_id: uniprot accession ID for protein B
        drug_to_fpt: dict maps drug ID to fingerprints
        protein_to_ligands: dict maps protein ID to set of drug IDs that bind that protein
        all_drug_ids: list of all drug IDs in drugs.csv
        num_iter: number of bootstrap iterations to perform
        seed: random seed 
    Returns:
        bootstrap p-value:
            p_bootstrap = sum(1(bootstrap Tsummary >= real Tsummary))) / num_iter

    """
    random.seed(seed)
    ligandsA = protein_to_ligands.get(proteinA_id)
    ligandsB = protein_to_ligands.get(proteinB_id)
    real_Tsummary = compute_tsummary(ligandsA, ligandsB, drug_to_fpt)

    num_greater_eq = 0  
    for i in range(num_iter):
        #sample w replacement
        #initially did this and was getting very low p-vals, realized this violates bootstrapping as a whole
        # random_ligands_A = set(random.choice(all_drug_ids) for i in range(len(ligandsA)))
        # random_ligands_B = set(random.choice(all_drug_ids) for i in range(len(ligandsA)))
        random_ligands_A = [random.choice(all_drug_ids) for i in range(len(ligandsA))]
        random_ligands_B = [random.choice(all_drug_ids) for i in range(len(ligandsB))]
        boot_Tsummary = compute_tsummary(random_ligands_A,random_ligands_B,drug_to_fpt)

        if boot_Tsummary>=real_Tsummary:
            num_greater_eq +=1
    return num_greater_eq / num_iter


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-n",type=int,help="num bootstrap iter.",default=500)
    parser.add_argument("-r",type=int, help="seed",default=214)
    parser.add_argument("drugs_file")
    parser.add_argument("targets_file")
    parser.add_argument("proteinA")
    parser.add_argument("proteinB")
    args = parser.parse_args()

    #load data
    drug_to_fpt = load_drugs_fpts(args.drugs_file)
    protein_to_ligands = load_protein_ligands(args.targets_file)
    all_drug_ids = list(drug_to_fpt)

    #compute p-val
    p_value = bootstrap_pvalue(args.proteinA,args.proteinB,drug_to_fpt,protein_to_ligands,all_drug_ids, num_iter= args.n,seed =args.r)
    print(p_value)