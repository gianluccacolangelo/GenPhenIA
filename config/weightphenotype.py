"""
En este script vamos a crear un .json donde cada fenotipo tiene un peso que es

n° genes asociados / n° genes totales en hpo
"""

##{{{ IMPORTACIONES
import json
import numpy as np
import networkx as nx
PATH = '/home/brainy/Desktop/1ercuatri2023/Tesis/GenPhenIA/'
import obonet
import pandas as pd
## }}}



## {{{
with open(f"{PATH}data/simulated/phenotype_gene_dict.json") as f:
    phenotype_gene_dict = json.load(f)

promiscuity_dict = {}

total_unique_genes = len(set(gene for genes in phenotype_gene_dict.values() for gene in genes))

for hpo_term, genes in phenotype_gene_dict.items():
    promiscuity = len(genes) / total_unique_genes
    promiscuity_dict[hpo_term] = promiscuity

# with open(f'{PATH}config/phen_promiscuity_dict.json', 'w') as f:
    # json.dump(promiscuity_dict, f, indent=4)

## }}}


## {{{
graph = obonet.read_obo(f'{PATH}data/hp.obo')
## }}}

## {{{

# Create a directed graph
G = nx.DiGraph()

# Add nodes and edges to the graph
for node in graph.nodes:
    G.add_node(node)
    for parent in graph.predecessors(node):
        G.add_edge(node, parent)

for node in G.nodes():
    if node not in promiscuity_dict:
        promiscuity_dict[node] = 'NaN'

# Calculate number of ancestors and children for each node
num_ancestors = {node: len(nx.ancestors(G, node)) for node in G.nodes()}
num_children = {node: len(nx.descendants(G, node)) for node in G.nodes()}
weights = {node: num_ancestors[node] / (num_ancestors[node] + num_children[node]) if num_ancestors[node] + num_children[node] > 0 else 0 for node in G.nodes()}


# Convert dictionaries to a pandas DataFrame
df = pd.DataFrame({'num_ancestors': pd.Series(num_ancestors),
                   'num_children': pd.Series(num_children),
                   'a/(a+c)': pd.Series(weights),
                   'frecuencia relativa': pd.Series(promiscuity_dict)})


# Save the DataFrame to a CSV file
df.to_csv('phen_properties.csv',sep='\t')



## }}}
