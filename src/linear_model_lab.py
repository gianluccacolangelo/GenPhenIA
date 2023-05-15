"""
En este script voy a ir definiendo las funciones que me permiten calcular y
evaluar los modelos
"""

## {{{ IMPORTACIONES
import numpy as np
import matplotlib.pyplot as plt
import phen_gen_weight_functions as pgw
import pandas as pd
PATH = "/home/brainy/Desktop/1ercuatri2023/Tesis/GenPhenIA/"
## }}}


## {{{
fen_observados = pgw.fen_observados_con_ruido(16,mph=0.3,iph=0.2,type_of_noise='normal')



def union_de_genes(set_of_phens):
    """
Esta función toma un conjunto de fenotipos y devuelvue la unión de todos los
genes que los causan.

Para eso recibe una lista de fenotipos de fen_observados_con_ruido.
Y para cada fenotipo scrappea en phenotypes_to_genes.txt para obtener los genes
que lo causan.

    """

    gene_ids = []
    df = pd.read_csv(f"{PATH}/data/phenotype_to_genes.txt",delimiter="\t")
    for hpo in set_of_phens:
        selected_rows = df[df["hpo_id"]==hpo]
        gene_ids.extend(selected_rows['ncbi_gene_id'].tolist())
    return set(gene_ids)


def calculate_gene_parameters(set_of_phens):

    genes = union_de_genes(set_of_phens)
    data = []
    i = 0

    for gene in genes:
        real_gene_phens = pgw.fen_reales_del_gen(gene)
        # calculate the parameters for this gene
        especificidad = pgw.especificidad_del_gen(set_of_phens,real_gene_phens)
        capitalidad = pgw.capitalidad_del_gen(set_of_phens,real_gene_phens)
        similaridad = pgw.similaridad_del_gen(set_of_phens,real_gene_phens)

        # add the gene and its parameters to the list
        data.append({'gene': gene,
                     'especificidad': especificidad,
                     'capitalidad': capitalidad,
                     'similaridad': similaridad})
        print(f"Calculando {i/len(genes)*100:.2f}%",end="\r")
        i+=1

    df = pd.DataFrame(data)

    df['total'] = df['especificidad'] + df['capitalidad'] + df['similaridad']

    # ordenamos el dataframe en orden descendente por el total
    df = df.sort_values('total', ascending=False)
    return df
## }}}
