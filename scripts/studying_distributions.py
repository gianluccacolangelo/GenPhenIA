"""
En este scritp vamos a estudiar
Especificidad(N)
Capitalidad(N)
Similaridad(N)

Donde N es el número máximo de fenotipos observados. Y lo vamos a hacer bajo
diferentes tipos de ruidos

La idea es ver cómo estos parámetros cambian a medida que aumentamos N para el
gen real. El gráfico final es un boxplot donde para cada N tenemos muchas
especificidades para cada gen real.
"""

##{{{ IMPORTACIONES
import numpy as np
import matplotlib.pyplot as plt
import random
import sys
sys.path.insert(0,'/home/brainy/Desktop/1ercuatri2023/Tesis/GenPhenIA/src')
import phen_gen_weight_functions as pgw
import linear_model_lab as lml
import scienceplots
import pandas as pd
##}}}


## {{{

list_of_genees = pgw.lista_de_genes()

def calculate_gene_parameters(mph,iph,type_of_noise,alpha,betha,gamma,gen_sample,fen_sample):
    """
Esta función calcula la especificidad, capitalidad y similaridad de gene
reales, dependiendo de cuánto pongamos en gen sample.
    """
    list_of_gens = np.random.choice(list_of_genees,gen_sample,replace=False)
    data = []
    noised_set = pgw.whats_your_set(mph,iph,type_of_noise)
    i=0
    for gene in list_of_gens:
        set_of_phens = pgw.fen_observados_con_ruido(gene,noised_set,fen_sample)
        real_gene_phens = pgw.fen_reales_del_gen(gene)
        # calculate the parameters for this gene
        especificidad = pgw.especificidad_del_gen(set_of_phens,real_gene_phens)
        capitalidad = pgw.capitalidad_del_gen(set_of_phens,real_gene_phens)
        similaridad = pgw.similaridad_del_gen(set_of_phens,real_gene_phens)

        # add the gene and its parameters to the list
        data.append({'gene': gene,
                     'especificidad':alpha*especificidad,
                     'capitalidad': betha*capitalidad,
                     'similaridad': gamma*similaridad,
                     'total':(especificidad+capitalidad+similaridad)})
        # print(f"Calculando {i/len(gene)*100:.1f}%",end="\r")
        # i+=1

    df = pd.DataFrame(data)

    # ordenamos el dataframe en orden descendente por el total
    df = df.sort_values('total', ascending=False)
    # Reseteamos el índice
    df = df.reset_index(drop=True)
    return df



def distributions(type_of_metric,
        mph,
        iph,
        type_of_noise,
        fen_sample,
        gen_sample,
        alpha,betha,gamma,
        list_of_genes=list_of_genees):

    df = calculate_gene_parameters(mph,iph,
            type_of_noise,
            alpha,betha,gamma,
            gen_sample,fen_sample)

    metric_dist = list(df[f'{type_of_metric}'])

    return metric_dist
## }}}


## {{{

especificidad_n = np.array([distributions('similaridad',.5,.5,'constant',
    fen_sample,500,1,1,1) for
        fen_sample in range(1,11)])

with plt.style.context(['science','ieee','nature']):
    fig,ax = plt.subplots()
    ax.boxplot(np.transpose(especificidad_n),meanline=True,
            flierprops=dict(markerfacecolor='r', markersize=0.4))
    ax.set_xlabel('Fenotipos observados totales',fontsize=4)
    ax.set_ylabel('Similaridad',fontsize=4)
    ax.set_title(f'500 HC para cada N. mph,iph=0.5',fontsize=4)
    plt.subplots_adjust(top=0.9,bottom=0.15,
            left=0.176,right=0.952,
            hspace=0.,wspace=0.1)

    ax.tick_params(axis='both', which='major', labelsize=4)
## }}}



