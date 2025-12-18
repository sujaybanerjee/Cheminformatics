import sys
import argparse
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

def load_node_info(protein_nodes_file):
    """
    Reads protein_nodes.csv data and builds mapping to uniprot_id and indication color
    Input: protein_nodes.csv: uniprot_accession,uniprot_id, indications
    Returns:
        id_map: dict mapping uniprot_accession to uniprot_id
        color_map: dict mapping uniprot_accession to color
    """
    df = pd.read_csv(protein_nodes_file)
    id_map = dict(zip(df["uniprot_accession"], df["uniprot_id"]))

    #specified color code
    colorcode = {
        "bp":"red",
        "bp;cholesterol":"green",
        "bp;cholesterol;diabetes":"blue",
        "bp;diabetes":"purple"}
    color_map = {}
    for i, row in df.iterrows():
        # color_map[row["uniprot_accession"]] = colorcode.get(row["indications"])
        color_map[row["uniprot_accession"]] = colorcode[row["indications"]]
    return id_map, color_map


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("edge_file")
    parser.add_argument("protein_nodes_file")
    parser.add_argument("output_png")
    args = parser.parse_args()

    #load edge list
    G = nx.read_edgelist(args.edge_file)

    #Relabeling and colorcoding
    id_map, color_map = load_node_info(args.protein_nodes_file)
    for accession in G.nodes():
        G.nodes[accession]["color"] = color_map.get(accession)
    G = nx.relabel_nodes(G, id_map)
    node_colors = [G.nodes[n]["color"] for n in G.nodes()]

    #draw graph
    plt.figure(figsize=(8,8))
    pos = nx.spring_layout(G,seed= 10)
    nx.draw_networkx(G,pos, node_color= node_colors,with_labels= True,font_size =8)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(args.output_png, dpi=150)
    plt.close()