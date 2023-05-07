"""
Este modulo almacena las funciones que le dan métrica a fenotipos y genotipos,
tanto en general, como para los pacientes en particular.
"""

## {{{ IMPORTACIONES
import numpy as np
## }}}


## {{{ functions

def raw_phen_to_genes(pacient_phen_list,phenotype_to_genes):
    """
Esta función toma una lista de fenotipos y devuelve una lista de todos los genes
candidatos. No pesa los fenotipos, no pesa los genes. Simplemente devuelve la
unión de los genes asociados a los fenotipos dados.
Para eso recorre la base de datos phenotype_to_genes.
    """

def pacient_gene_score(pacient_genes_list):
    """
Esta función toma la lista de raw_phen_to_genes del paciente y le asigna un score
a los elementos en base a su capitalidad y su especificidad con respecto a
pacient_phen_list
    """

def pacient_phen_score(pacient_phen_list):
    """
Esta función toma la lista de fenotipos y devuelve un diccionario (u otra
clase) con sus scores.
    """

def general_gene_score(genes_to_phenotype):
    """
Esta función toma el archivo genes_to_phenotype.txt y ...
acá hay que hacer redes y centralidades
    """

def general_phen_score(phenotype_to_genes):
    """
Esta función toma el archivo phenoype_to_genes.txt y le asigna un score a todos
los fenotipos según la relación 1/n° genes que lo causan. De modo que aquél
fenotipo que tiene una relación 1:1 tiene el máximo score.
    """



## }}}
