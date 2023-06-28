"""
En este script vamos a construir y analizar una red gen-fenotipo utilizando datos de Hu-
man Phenotype Ontology (HPO) y las bases de datos de OMIM, con el objetivo
de identificar comunidades de genes asociados con fenotipos similares
"""


## {{{
import networkx as nx
import json
import community as community_louvain
import numpy as np
import pandas as pd
import sys
sys.path.insert(0,'/home/brainy/Desktop/1ercuatri2023/Tesis/GenPhenIA/src')
import phen_gen_weight_functions as pgw
from networkx.algorithms import bipartite
import matplotlib.pyplot as plt
import community as community_louvain


PATH = "/home/brainy/Desktop/1ercuatri2023/Tesis/GenPhenIA/"
## }}}


## {{{
with open(f"{PATH}data/genes_to_phenotype.txt",'r') as f:
    genes_to_phenotypes = pd.read_csv(f, sep='\t', header=None)

all_genes = list(genes_to_phenotypes[1].unique())

## }}}


## {{{

#with the function pgw.gene_diseases(gene_symbol) create a dictionary of genes
#and their diseases, noting that a gene can cause multiple diseases.

genes_diseases = {}
i=0
for gene in all_genes:
    gene_diseases = pgw.gene_diseases(gene)
    for disease in gene_diseases:
       classification = pgw.disease_classification(disease,1)
       if classification is not None:
            if gene not in genes_diseases:
                genes_diseases[gene] = set()
            genes_diseases[gene].add(classification[1])

## }}}


## {{{
# Convert sets back to lists
for gene in genes_diseases:
    genes_diseases[gene] = list(genes_diseases[gene])
#Now create a bipartite network of genes and diseases from the dictionary created

B = nx.Graph()

B.add_nodes_from(genes_diseases.keys(), bipartite=0) # add genes as one partition

all_diseases = set([disease for diseases in genes_diseases.values() for disease in diseases])

B.add_nodes_from(all_diseases, bipartite=1) # add diseases as second partition


for gene, diseases in genes_diseases.items():
    for disease in diseases:
        B.add_edge(gene, disease)


# assert bipartite.is_bipartite(B)

G_genes = bipartite.projected_graph(B, genes_diseases.keys())
G_diseases = bipartite.projected_graph(B, all_diseases)
## }}}


## {{{

# Create weighted projected graphs
G_genes = bipartite.weighted_projected_graph(B, genes_diseases.keys())
G_diseases = bipartite.weighted_projected_graph(B, all_diseases)
## }}}


## {{{
partition = community_louvain.best_partition(G_diseases)


## }}}


## {{{
threshold = 3

# Iterate through the edges of the graph
for u, v, data in list(G_genes.edges(data=True)):
    # If the weight is below the threshold, remove the edge
    if data['weight'] < threshold:
        G_genes.remove_edge(u, v)

# Similarly for the diseases graph
for u, v, data in list(G_diseases.edges(data=True)):
    if data['weight'] < threshold:
        G_diseases.remove_edge(u, v)

## }}}



## {{{
# Get the degree of the nodes in the bipartite graph B
degree_dict = dict(B.degree(all_diseases))
degrees = [degree_dict[node] for node in G_diseases.nodes()]

# draw the graph
pos = nx.spring_layout(G_diseases)
plt.figure(figsize=(18, 18))
# plt.axis('off')

labels = {node: pgw.disease_classification(node)[0] for node in G_diseases.nodes()}

for key in labels.keys():
    labels[key] = labels[key].replace('Rare', '').replace('genetic',
            '').replace(' ','\n').strip()


# Draw nodes in partition-specific colors

# Draw nodes in partition-specific colors
for community in set(partition.values()):
    node_list = [node for node in partition.keys() if partition[node] == community]
    node_sizes = [degree_dict[node] * 10 for node in node_list]  # Get degrees for the current community nodes
    nx.draw_networkx_nodes(G_diseases, pos, node_list,
                           node_color=plt.cm.nipy_spectral(community / len(set(partition.values()))),
                           node_size=node_sizes,alpha=.5)  # Use the new node_sizes list
nx.draw_networkx_edges(G_diseases, pos, alpha=0.15)
label_sizes = {node: 10 + np.log(degree_dict[node])  for node in G_diseases.nodes()}  # Create a dictionary with label sizes proportional to node degree (adjust the divisor for better results)

# nx.draw_networkx_labels(G_diseases, pos,
        # labels,font_weight='bold',font_size=label_sizes)
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman']})
plt.rc('text', usetex=True)
plt.rcParams['font.weight'] = 'bold'



for node, (x, y) in pos.items():
    plt.text(x, y,
             s=labels[node],
             fontsize=label_sizes[node],
             ha='center', va='center',
             fontweight='bold')
plt.show()

#In this code, community_louvain.best_partition(G) is used to find the community structure of the graph G that maximizes the modularity, using the Louvain method. Each community will be colored differently on the graph. The labels of the nodes are the gene names. The position of each node (gene) is determined using the spring layout, which treats edges as springs holding nodes close, while treating nodes as repelling objects. Adjusting the spring layout can help make the graph more understandable.






## }}}



## {{{
pos = nx.bipartite_layout(B, genes_diseases.keys())
plt.figure(figsize=(10, 10))
nx.draw(B, pos, with_labels=True,
        nodelist=genes_diseases.keys(),
        node_color='blue',
        # node_size=[B.degree(gene) * 500 for gene in genes_diseases.keys()],
        alpha=0.2)

nx.draw_networkx_nodes(B, pos,
                       nodelist=all_diseases,
                       node_color='red',
                       node_size=[B.degree(disease) * 50 for disease in all_diseases],
                       alpha=0.2)


nx.draw_networkx_edges(B, pos, alpha=0.01)

nx.draw_networkx_labels(B, pos)
plt.show()
## }}}
