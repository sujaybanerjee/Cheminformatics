README.md

Chemoinformatics: Protein Similarity via Ligand Chemistry

This project applies chemoinformatics methods to compare proteins based on the chemical similarity of their ligands rather than sequence or structure alone. Using DrugBank data, ligand fingerprints, and Tanimoto similarity, the pipeline identifies statistically significant relationships between proteins and visualizes them as networks.

Motivation

Proteins that are dissimilar in sequence or structure may still bind similar ligands, while closely related proteins may not share ligands at all. By comparing proteins through the chemistry of their binding partners, this project follows a ligand-centric approach inspired by the Similarity Ensemble Approach (SEA).

Methods Overview
	•	Chemical representation
	•	Parsed SMILES strings and molecular fingerprints from DrugBank
	•	Compared ligands using the Tanimoto (Jaccard) coefficient
	•	Ligand similarity analysis
	•	Computed all pairwise Tanimoto scores between drugs
	•	Evaluated whether chemically similar ligands share protein targets
	•	Protein similarity via ligand sets
	•	Defined a ligand-set similarity score by summing Tanimoto values above a threshold (Tc > 0.5)
	•	Assessed significance using bootstrap resampling to generate empirical p-values
	•	Network construction
	•	Built a protein–protein network where edges represent statistically significant ligand-set similarity
	•	Visualized the network with NetworkX, coloring nodes by disease indication

Key Scripts
	•	tanimoto.py
Computes pairwise Tanimoto similarities between drugs and labels shared targets.
	•	pvalue.py
Calculates bootstrap p-values for ligand-set similarity between two proteins.
	•	networkgen.py
Generates a protein interaction edgelist based on significant bootstrap results.
	•	plot_graph.py
Visualizes the protein similarity network using NetworkX and Matplotlib.

Outputs
	•	Histograms of ligand similarity distributions
	•	Bootstrap p-values for protein pairs
	•	Protein similarity network highlighting shared pharmacology across disease areas

Technologies Used
	•	Python
	•	NetworkX
	•	NumPy / Pandas
	•	Matplotlib
	•	DrugBank data
