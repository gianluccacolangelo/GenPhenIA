"""
Este modulo almacena las funciones que le dan peso a fenotipos y genotipos,
tanto en general, como para los pacientes en particular.
"""

## {{{ IMPORTACIONES
import numpy as np
import json
PATH = '/home/brainy/Desktop/1ercuatri2023/Tesis/GenPhenIA/'
## }}}


## {{{ functions

def raw_phen_to_genes(pacient_phen_list,phenotype_to_genes):
    """
Esta función toma una lista de fenotipos y devuelve una lista de todos los genes
candidatos. No pesa los fenotipos, no pesa los genes. Simplemente devuelve la
unión de los genes asociados a los fenotipos dados.
Para eso recorre la base de datos phenotype_to_genes.
    """

## {{{ funciones de peso de genes candidatos
def especificidad_del_gen(fenotipos_observados,fenotipos_del_gen):
    """
fenotipos_observados y fenotipos_del_gen tienen que ser sets (conjuntos), esta
función devolverá la fracción de fenotiops_del_gen que están en fenotipos_observados

Recordar que siempre hablamos del gen candidato para el conj de fenotipos obs.
    """
    return len(fenotipos_observados.intersection(
        fenotipos_del_gen))/len(fenotipos_del_gen)


def capitalidad_del_gen(fenotipos_observados,fenotipos_del_gen):
    """
esta función fevolverá la fracción de fenotiops_observados que están en fenotipos_del_gen

Recordar que siempre hablamos del gen candidato para el conj de fenotipos obs.
    """
    return len(fenotipos_observados.intersection(
        fenotipos_del_gen))/len(fenotipos_observados)


def similaridad_del_gen(fenotipos_observados,fenotipos_del_gen):
    """
Esta función devolverá la intersección de fenotipos observados y fenotipos del
gen sobre la unión de ambos.

Recordar que siempre hablamos del gen candidato para el conj de fenotipos obs.
    """
    return len(fenotipos_observados.intersection(
        fenotipos_del_gen))/len(fenotipos_observados.union(fenotipos_del_gen))


## }}}


## {{{ funciones que toman los conj de fen. observados y fen. reales
"""
Hay dos formas de hacer esto, una es la propuesta abajo, que es tomar los
fenotipos "observados" con los sets simulados con ruido

La otra es no usar los sets ruidosos y usar el archivo
phenotype_to_genes, que tiene la ventaja de que es más amplio, en el sentido
que hay links entre fenotipos muuuy generales y los genes, que es algo que
suele suceder entre los médicos, que a veces te tiran "anomalía de
crecimiento", y puede ser básicamente cualq cosa.

La tercera opción sería combinar ambas cosas. Pero paso a paso.
"""


def fen_reales_del_gen(gen_id,
        db=f'{PATH}data/simulated/gene_phenotype_dict.json'):
    # gene_id = str(gen_id)
    with open(db,'r') as file:
        gene_phenotype_dict = json.load(file)
    fen_reales = gene_phenotype_dict[str(gen_id)]
    return set(fen_reales)

def fen_observados_con_ruido(gen_id,
        mph=0.1,
        iph=0.1):
    """
mph e iph corresponden a los dos parámetros de ruido: missing phenotypes e
incorrect phenotypes, respectivamente, el que elijamos ahí va a ser tomado de
la base de datos.
    """
    db=f'{PATH}data/simulated/simulated_set_mph{mph}_iph{iph}.json'
    with open(db,'r') as file:
        noised_gene_phenotype_dict = json.load(file)
    fen_observados = dict(noised_gene_phenotype_dict)[str(gen_id)]
    return set(fen_observados)


## }}}




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
