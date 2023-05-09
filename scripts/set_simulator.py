"""
Este script contiene las funciones que me van a permitir generar set simulados
a partir de la base de datos genes_to_phenotype y phenotype_to_genes.
"""

## {{{ Importaciones
import numpy as np
import random
## }}}


def txt_to_dict(filename):
    """
Esta función toma ya sea genes_to_phenotype.txt o phenotype_to_genes.txt y la
convierte en un diccionario con la siguiente forma
{'gene1':['phen1','phen2','phen3'], gene2: ....}
    """


def generate_incorrect_phenotypes(phenotypes,num_incorrect):
    """
Genera una lista de genes incorrectos a partir de una lista de fenotipos de un
dado gen

phenotypes (List[str]): Una lista de fenotipos asociados con un gen específico
num_incorrect (int): El número de fenotipos incorrectos a generar

devuelve: una lista de {num incorrect} fenotipos incorrectos
    """
    all_phenotypes = set()
    for phenlist in gene_phenotype_dict.values():
        for phen in phenlist:
            all_phenotypes.add(p)

    phenotypes_not_in_current_gene = []
    for phen in all_phenotypes:
        if phen not in phenotypes:
            phenotypes_not_in_current_gene.append(phen)

    incorrect_phenotypes = random.sample(phenotypes_not_in_current_gene,
            num_incorrect)
    return incorrect_phenotypes


