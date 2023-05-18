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

df_phen_to_gen = pd.read_csv(f"{PATH}/data/phenotype_to_genes.txt",delimiter="\t")
def union_de_genes(set_of_phens):
    """
Esta función toma un conjunto de fenotipos y devuelvue la unión de todos los
genes que los causan.

Para eso recibe una lista de fenotipos de fen_observados_con_ruido.
Y para cada fenotipo scrappea en phenotypes_to_genes.txt para obtener los genes
que lo causan.

    """
    selected_rows = df_phen_to_gen[df_phen_to_gen["hpo_id"].isin(set_of_phens)]
    gene_ids = selected_rows['ncbi_gene_id'].tolist()
    return set(gene_ids)

def calculate_gene_parameters(set_of_phens):

    genes = union_de_genes(set_of_phens)
    data = []
    i = 0

    for gene in genes:
        i+=1

        real_gene_phens = pgw.fen_reales_del_gen(gene)
        # calculate the parameters for this gene
        especificidad = pgw.especificidad_del_gen(set_of_phens,real_gene_phens)
        capitalidad = pgw.capitalidad_del_gen(set_of_phens,real_gene_phens)
        similaridad = pgw.similaridad_del_gen(set_of_phens,real_gene_phens)

        # add the gene and its parameters to the list
        data.append({'gene': gene,
                     'especificidad': especificidad,
                     'capitalidad': capitalidad,
                     'similaridad': similaridad,
                     'total':(especificidad+capitalidad+similaridad)/3})
        print(f"Calculando {i/len(genes)*100:.1f}%",end="\r")

    df = pd.DataFrame(data)

    # ordenamos el dataframe en orden descendente por el total
    df = df.sort_values('total', ascending=False)
    # Reseteamos el índice
    df = df.reset_index(drop=True)
    return df
## }}}


## {{{ Evaluando para todos los genes

def model_evaluating(mph,iph,type_of_noise):
    """
A esta función le damos un mph y un iph y evalúa el modelo para ese set
simulado, calculando el ranking en la importancia para el gen real dado un
conjunto de fenotipos.
    """
    list_of_genes = pgw.lista_de_genes()
    #metrics va a almacenar la posición en el ranking de importancia del gen
    # real en nuestro modelo
    metrics = []
    i=1
    noised_set = pgw.whats_your_set(mph,iph,"normal")
    for gene in list_of_genes:
        print(f"\nCalculando para gen {i} de {len(list_of_genes)}  ")
        #Los fenotipos observados con ruido para un dado gen
        fen_observados = pgw.fen_observados_con_ruido(gene,
                noised_set)

        #Calculamos los parámetros esp. cap. y sim. para la unión de genes
        # posibles que causan esos fenotipos
        df = calculate_gene_parameters(fen_observados)

        #Esto agrega a metrics la posición en el índice rankeado del gen real
        # entre los miles posibles
        try:
            metrics.append(df.loc[df['gene']==int(gene)].index[0])
            print(f"                      ranking = {metrics[i-1]+1}\n")
            i+=1
        except:
            continue


    return metrics



## }}}


