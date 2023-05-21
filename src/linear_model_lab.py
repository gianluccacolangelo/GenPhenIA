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

def calculate_gene_parameters(set_of_phens,alpha,betha,gamma):
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

        # add the gene and its parameters to the list
        data.append({'gene': gene,
                     'especificidad':especificidad,
                     'capitalidad': capitalidad,
                     'similaridad': similaridad,
                     'total':(alpha*especificidad+betha*capitalidad+gamma*similaridad)})
        print(f"Calculando {i/len(genes)*100:.1f}%",end="\r")

    df = pd.DataFrame(data)

    # ordenamos el dataframe en orden descendente por el total
    df = df.sort_values('total', ascending=False)
    # Reseteamos el índice
    df = df.reset_index(drop=True)
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
        df = calculate_gene_parameters(fen_observados,alpha,betha,gamma)

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





## {{{
# alpha,betha,gamma = (0,1,0)

# for alpha,betha,gamma in [(1,0,0)]:
    # results = {}
# #alpha -> especificidad
# #betha -> capitalidad
# #gamma -> similaridad
# # for type_of_noise in ['normal','constant','random','gold_standard']:
    # for type_of_noise in ['constant','gold_standard']:
        # if type_of_noise == 'gold_standard':
            # mph_iph_metrics = []
            # for fen_samples in range(10):
                # list_of_gen_rankings = model_evaluating(0.1,0.1,type_of_noise,
                            # fen_samples+1,
                            # 500,
                            # alpha,
                            # betha,
                            # gamma)
                # mph_iph_metrics.append(percent_below_x(list_of_gen_rankings,10))
            # results[f"clean_set"] = mph_iph_metrics


        # elif type_of_noise == 'normal' or type_of_noise == 'constant':
            # for mph in [0.1]:
                # for iph in [0.1]:
                    # mph_iph_metrics = []
                    # for fen_samples in range(10):
                        # list_of_gen_rankings = model_evaluating(mph,iph,
                                    # type_of_noise,
                                    # fen_samples+1,
                                    # 500,
                                    # alpha,
                                    # betha,
                                    # gamma)
                        # mph_iph_metrics.append(percent_below_x(list_of_gen_rankings,10))
                    # results[f"{type_of_noise}: mph={mph}, iph={iph}"] = mph_iph_metrics

        # elif type_of_noise == 'random':
            # mph_iph_metrics = []
            # for fen_samples in range(15):
                # list_of_gen_rankings = model_evaluating(0.1,0.1,type_of_noise,fen_samples+1,500)
                # mph_iph_metrics.append(percent_below_x(list_of_gen_rankings,10))
            # results[f"random"] = mph_iph_metrics


        # top_10_metrics = pd.DataFrame(results)

        # with open(f"{PATH}output/top_10_metrics_{alpha}_{betha}_{gamma}.csv", "w") as f:
            # top_10_metrics.to_csv(f)

## }}}


## {{{
# with open(f"{PATH}output/top_10_metrics_1_0_0.csv", "r") as f:
    # top_10_metrics = pd.read_csv(f)
## }}}


## {{{
# x = np.linspace(0, 10, 100)

# # Create a set of line styles to use (you can use any suitable line styles)
# styles = ['-', '--', '-.', ':']

# # Create a colormap
# cmap = plt.get_cmap('tab10')

# with plt.style.context(['science','ieee','nature']):
    # fig, ax1 = plt.subplots()
    # plt.rcParams["text.usetex"] = True
    # for i, label in enumerate(top_10_metrics.columns[1:2]):
        # ax1.plot(np.arange(1,11), top_10_metrics[label], label=label, color=cmap(i), linestyle=styles[i % len(styles)])
    # ax1.plot(np.arange(1,11), top_10_metrics['clean_set'], label='clean_set',
            # color='black')
    # ax1.legend(loc='lower right', fontsize=4)
    # ax1.set_xlabel('Total observed phenotypes (N)', fontsize=4)
    # ax1.set_ylabel('Accuracy',fontsize=4)
    # ax1.set_title('Accuracy of the model for $\\alpha=1, \\beta=\\gamma=0$', fontsize=4)
    # ax1.tick_params(axis='both', which='major', labelsize=3)
    # # ax1.axhline(y=0.9, color='green', linestyle='--', linewidth=1,alpha=0.5)
    # # ax1.axvline(x=5, color='green', linestyle='--', linewidth=1,alpha=0.5)
    # ax1.set_ylim(0,1)
    # plt.subplots_adjust(top=0.9,bottom=0.15,
            # left=0.176,right=0.952,
            # hspace=0.,wspace=0.1)
    # for i, label in enumerate(top_10_metrics.columns[1:5]):
        # ax2.plot(np.arange(0,15), top_10_metrics[label], label=label, color=cmap(i), linestyle=styles[i % len(styles)])
    # ax2.plot(np.arange(0,15), top_10_metrics['clean_set'], label='clean_set',
            # color='black')
    # ax2.legend(loc='lower right', fontsize=2)
    # ax2.set_xlabel('Total observed phenotypes', fontsize=4)
    # ax2.set_ylabel('Accuracy',fontsize=4)
    # ax2.set_title('Accuracy of the model for different noise levels', fontsize=4)
    # ax2.tick_params(axis='both', which='major', labelsize=4)
    # ax2.axhline(y=0.9, color='green', linestyle='--', linewidth=1,alpha=0.5)
    # ax2.axvline(x=5, color='green', linestyle='--', linewidth=1,alpha=0.5)

# plt.show()
## }}}


##{{{

# best_performance = -np.inf
# best_weights = None

# for alpha in np.arange(0,1,0.1):
    # for betha in np.arange(0,1,0.1):
        # for gamma in np.arange(0,1,0.1):
            # print(f"alpha={alpha}, betha={betha}, gamma={gamma}")
            # performance = percent_below_x(model_evaluating(0.2,0.2,'constant',3,250,alpha,betha,gamma),10)
            # if performance > best_performance:
                # best_performance = performance
                # best_weights = (alpha,betha,gamma)
## }}}



