"""
En este script voy a ir definiendo las funciones que me permiten calcular y
evaluar los modelos
"""

## {{{ IMPORTACIONES
import numpy as np
import matplotlib.pyplot as plt
import phen_gen_weight_functions as pgw
import pandas as pd
import random
import scienceplots
PATH = "/home/brainy/Desktop/1ercuatri2023/Tesis/GenPhenIA/"
## }}}


## {{{

df_phen_to_gen = pd.read_csv(f"{PATH}/data/phenotype_to_genes.txt",delimiter="\t")
phen_to_genes = df_phen_to_gen.groupby('hpo_id')['ncbi_gene_id'].apply(set).to_dict()

def union_de_genes(set_of_phens):
    """
Esta función toma un conjunto de fenotipos y devuelvue la unión de todos los
genes que los causan.

Para eso recibe una lista de fenotipos de fen_observados_con_ruido.
Y para cada fenotipo scrappea en phenotypes_to_genes.txt para obtener los genes
que lo causan.

    """
    gene_ids = set.union(*(phen_to_genes.get(phen, set()) for phen in set_of_phens))
    return gene_ids

def calculate_gene_parameters(set_of_phens,alpha,betha,gamma,n_metrica=1,nueva_metrica=False):
    """
Esta función toma un conjunto de fenotipos observados, y calcula:
    especificidad
    capitalidad
    similaridad

    para cada uno de los genes candidatos y devuelve una lista ordenada por
    aquellos genes candidatos que más suman esas métricas pesadas por alpha
    betha y gamma.
    """

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
        nueva_metrica = pgw.parametro(set_of_phens,real_gene_phens,n_metrica)

        # add the gene and its parameters to the list
        data.append({'gene': gene,
                     'especificidad':especificidad,
                     'capitalidad': capitalidad,
                     'similaridad': similaridad,
                     'nueva_metrica':nueva_metrica,
                     'total':(alpha*especificidad+betha*capitalidad+gamma*similaridad)})
        print(f"Calculando {i/len(genes)*100:.1f}%",end="\r")

    df = pd.DataFrame(data)

    # ordenamos el dataframe en orden descendente por el total
    if nueva_metrica == False:
        df = df.sort_values('total', ascending=False)
        # Reseteamos el índice
        df = df.reset_index(drop=True)
        print("sorting by total")
    else:
        df = df.sort_values('nueva_metrica', ascending=False)
        df = df.reset_index(drop=True)
        print("sorting by new metric")
    return df
## }}}


## {{{ Evaluando para todos los genes

list_of_gens = pgw.lista_de_genes()

def model_evaluating(mph,iph,
        type_of_noise,
        fen_sample,
        gen_sample,
        alpha,
        betha,
        gamma,
        nueva_metrica=False,
        n_metrica=1,
        list_of_genes=list_of_gens):
    """
A esta función le damos un mph, iph, y tipo de ruido y evalúa el modelo para ese
set simulado, para cada uno de los genes reales, de una seleccion gen_sample al
azar de la lista total de genes.

fen_sample (N) es el número máximo de fenotipos observados que vamos a obtener.

finalmente devuelve una lista con la posición en el ranking de cada uno de los
genes reales evaluados
    """
    #metrics va a almacenar la posición en el ranking de importancia del gen
    # real en nuestro modelo
    list_of_genes = np.random.choice(list_of_genes,gen_sample,replace=False)
    metrics = []
    i=1
    noised_set = pgw.whats_your_set(mph,iph,type_of_noise)
    for gene in list_of_genes:
        print(f"\nCalculando para gen {i} de {len(list_of_genes)}  ")
        #Los fenotipos observados con ruido para un dado gen
        fen_observados = pgw.fen_observados_con_ruido(gene,
                noised_set,
                fen_sample)

        #Calculamos los parámetros esp. cap. y sim. para la unión de genes
        # posibles que causan esos fenotipos

        df = calculate_gene_parameters(fen_observados,
                alpha,betha,gamma,
                nueva_metrica=nueva_metrica,
                n_metrica=n_metrica)

        #Esto agrega a metrics la posición en el índice rankeado del gen real
        # entre los miles posibles
        try:
            metrics.append(df.loc[df['gene']==int(gene)].index[0]+1)
            print(f"                      ranking = {metrics[i-1]}\n")
            i+=1
        except:
            continue


    return metrics


def percent_below_x(lst,x):
    count = sum(1 for i in lst if i <= x)
    return (count / len(lst))
## }}}

