import pandas as pd

def load_drugs_fpts(drugs_file):
    """
    Loads drug fingerprints from drugs.csv
    Inputs:
        drugs_file: Path to drugs.csv containing columns:
            drug_id
            fingerprint(maccs)
    Returns:
        dict: drug_id to set[int] of fingerprint keys
    """
    df = pd.read_csv(drugs_file)
    drug_to_fpt = {}
    for _, row in df.iterrows():
        drug_id = row['db_id']
        fpt_string = str(row['maccs']).strip()
        #if missing
        if fpt_string == "":
            drug_to_fpt[drug_id] = set()
        else:
            # fpt_set = set(int(x) for x in fpt_string.split())
            fpt_set = set(fpt_string.split())
            drug_to_fpt[drug_id] = fpt_set
    return drug_to_fpt



def load_protein_ligands(targets_file):
    """
    Loads protein to drug relationships from targets.csv
    Inputs:
        targets_file: path to targets.csv:
            db_id
            uniprot_accession
    Returns:
        dict: protein_id to set of db_ids that bind it
    """
    df = pd.read_csv(targets_file)
    protein_to_drugs = {}
    for _,row in df.iterrows():
        protein = row['uniprot_accession']
        drug = row['db_id']
        if protein not in protein_to_drugs:
            protein_to_drugs[protein] = set()
        protein_to_drugs[protein].add(drug)
    return protein_to_drugs


def tanimoto(fpt_a, fpt_b):
    """
    Computes Tanimoto similarity between fingerprint sets
    Inputs:
        fpt_a - fingerprint a
        fpt_b - fingerprint b
    Returns:
        float tanimoto score
    """
    #|fpt(mol_a)and fpt(mol_b)| / |fpt(mol_a) or fpt(mol_b)|
    intersection = len(fpt_a & fpt_b)
    union = len(fpt_a | fpt_b)
    if union == 0:
        return 0.0
    return intersection / union