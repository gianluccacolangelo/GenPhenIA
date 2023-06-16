"""
En este script vamos a crear un .json donde cada fenotipo tiene un peso que es

n° genes asociados / n° genes totales en hpo
"""

##{{{ IMPORTACIONES
import json
import numpy as np
PATH = '/home/brainy/Desktop/1ercuatri2023/Tesis/GenPhenIA/'
## }}}



## {{{
with open(f"{PATH}data/simulated/phenotype_gene_dict.json") as f:
    phenotype_gene_dict = json.load(f)

promiscuity_dict = {}

total_unique_genes = len(set(gene for genes in phenotype_gene_dict.values() for gene in genes))

for hpo_term, genes in phenotype_gene_dict.items():
    promiscuity = 1 -  len(genes) / total_unique_genes
    promiscuity_dict[hpo_term] = promiscuity

with open(f'{PATH}config/phen_promiscuity_dict.json', 'w') as f:
    json.dump(promiscuity_dict, f, indent=4)

## }}}
