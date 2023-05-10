"""
Este script contiene las funciones que me van a permitir generar set simulados
a partir de la base de datos genes_to_phenotype y phenotype_to_genes.
"""

## {{{ Importaciones
import numpy as np
import random
import csv
import json
## }}}


## {{{ txt_to_dict

#REEMPLAZAR POR EL PATH DONDE TENES TU REPOSITORIO GenPhenIA
PATH = '/home/brainy/Desktop/1ercuatri2023/Tesis/GenPhenIA/'

def txt_to_dict(filename):
    """
Esta función toma ya sea genes_to_phenotype.txt o phenotype_to_genes.txt y la
convierte en un diccionario con la siguiente forma
{'gene1':['phen1','phen2','phen3'], gene2: ....} almacenándolo en formato .json
    """
    gene_phenotype_dict = {}
    with open(filename,'r') as file:
        reader = csv.reader(file,delimiter='\t')
        next(reader)
        for row in reader:
            gene_id, gene_symbol, hpo_id, hpo_name = row
            if gene_id not in gene_phenotype_dict:
                gene_phenotype_dict[gene_id] = []
            gene_phenotype_dict[gene_id].append(hpo_id)

    with open(f'{PATH}/data/simulated/gene_phenotype_dict.json','w') as file:
        json.dump(gene_phenotype_dict, file)
## }}}

## {{{ opening dict
with open(f'{PATH}data/simulated/gene_phenotype_dict.json', 'r') as file:
    gene_phenotype_dict = json.load(file)
## }}}

## {{{ generate_incorrect_phenotypes
def generate_incorrect_phenotypes(phenotypes,num_incorrect,all_phenotypes):
    """
Genera una lista de genes incorrectos a partir de una lista de fenotipos de un
dado gen

phenotypes (List[str]): Una lista de fenotipos asociados con un gen específico
num_incorrect (int): El número de fenotipos incorrectos a generar

devuelve: una lista de {num incorrect} fenotipos incorrectos
    """
    phenotypes_not_in_current_gene = []
    for phen in all_phenotypes:
        if phen not in phenotypes:
            phenotypes_not_in_current_gene.append(phen)

    incorrect_phenotypes = random.sample(phenotypes_not_in_current_gene,
            num_incorrect)
    return incorrect_phenotypes
## }}}

## {{{ genphen_simulator

def genphen_simulator(missing_phens=0.1,
        incorrect_phens=0.1,
        n_samples_per_gene=10,
        genphen_db=f'{PATH}/data/simulated/genes_to_phenotype.txt'):
    """
Esta función toma la base de datos genes_to_phenotype.txt y genera un set
simulado a partir del porcentaje de missing_phns e incorrect_phens que le
demos,y genera n_samples por gen tomando de la db genes_to_phenotype.txt.
    """
    all_phenotypes = set()
    for phenlist in gene_phenotype_dict.values():
        for phen in phenlist:
            all_phenotypes.add(phen)

    simulated_data = []
    n=0
    for gene, phenotypes in gene_phenotype_dict.items():
        missing_num = int(np.round(len(phenotypes) * missing_phens))
        incorrect_num = int(np.round(len(phenotypes) * incorrect_phens))
        for i in range(n_samples_per_gene):

            # Generate incomplete phenotypes
            missing_phenotypes = random.sample(phenotypes, missing_num)

            # Generate incorrect phenotypes
            incorrect_phenotypes = generate_incorrect_phenotypes(phenotypes,
                    incorrect_num,
                    all_phenotypes)

            # Create modified phenotype list
            simulated_phenotypes = [p for p in phenotypes if p not in missing_phenotypes] + incorrect_phenotypes
            random.shuffle(simulated_phenotypes)

            # Add modified phenotype list and associated gene to the simulated dataset
            simulated_data.append((gene, simulated_phenotypes))
        print(f'Generando base: {n/4921*100}%')
        n+=1

    with open(f'{PATH}/data/simulated/simulated_set.json','w') as file:
        json.dump(simulated_data, file)
## }}}
