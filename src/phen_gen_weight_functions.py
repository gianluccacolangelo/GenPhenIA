"""
Este modulo almacena las funciones que le dan peso a fenotipos y genotipos,
tanto en general, como para los pacientes en particular.
"""

## {{{ IMPORTACIONES
import numpy as np
import json
import pandas as pd
import matplotlib.pyplot as plt
PATH = '/home/brainy/Desktop/1ercuatri2023/Tesis/GenPhenIA/'
## }}}




## {{{ funciones de peso de genes candidatos
def especificidad_del_gen(fenotipos_observados,fenotipos_del_gen):
    """
fenotipos_observados y fenotipos_del_gen tienen que ser sets (conjuntos), esta
función devolverá la fracción de fenotiops_del_gen que están en fenotipos_observados

Recordar que siempre hablamos del gen candidato para el conj de fenotipos obs.
    """
    j = fenotipos_observados.intersection(fenotipos_del_gen)
    i = set(fenotipos_observados)-j
    k = set(fenotipos_del_gen)-j
    # if type_of_func == "log":
        # result = len(j)/(np.log10(1+len(k))+len(j))
    # elif type_of_func == "exp":
        # result = len(j)/(10**(len(k))+len(j))
    # elif type_of_func == "lineal":
        # result = len(j)/(len(k)+len(j))
    return len(j)/(np.log10(len(k)+1)+len(j))


def capitalidad_del_gen(fenotipos_observados,fenotipos_del_gen):
    """
esta función fevolverá la fracción de fenotiops_observados que están en fenotipos_del_gen

Recordar que siempre hablamos del gen candidato para el conj de fenotipos obs.
    """
    j = fenotipos_observados.intersection(fenotipos_del_gen)
    i = len(set(fenotipos_observados)-j)
    k = len(set(fenotipos_del_gen)-j)
    j = len(j)
    return j/(i+j)

def similaridad_del_gen(fenotipos_observados,fenotipos_del_gen):
    """
Esta función devolverá la intersección de fenotipos observados y fenotipos del
gen sobre la unión de ambos.

Recordar que siempre hablamos del gen candidato para el conj de fenotipos obs.
    """
    j = fenotipos_observados.intersection(fenotipos_del_gen)
    i = len(set(fenotipos_observados)-j)
    k = len(set(fenotipos_del_gen)-j)
    j = len(j)
    return j/(j+i+np.log10(1+k))

def parametro(fenotipos_observados,fenotipos_del_gen,nuevo_parametro):
    """
Esta función está para probar nuevas métricas, donde nuevo_parametro es un
entero de la nueva métrica
    """
    j = fenotipos_observados.intersection(fenotipos_del_gen)
    i = set(fenotipos_observados)-j
    k = set(fenotipos_del_gen)-j
    if nuevo_parametro==1:
        result = len(j)-len(i)-len(k)
    elif nuevo_parametro==2:
        result = len(j)-len(i)
    elif nuevo_parametro==3:
        result = len(j)/(1+len(i)+len(k))
    elif nuevo_parametro==4:
        result = len(j)
    elif nuevo_parametro==5:
        result = -len(i)
    elif nuevo_parametro==6:
        result = -len(k)
    elif nuevo_parametro == 7:
        result = -len(i)-len(k)
    return result


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


# Abrimos gene_phenotype_dict
db=f'{PATH}data/simulated/gene_phenotype_dict.json'
with open(db,'r') as file:
    gene_phenotype_dict = json.load(file)

def fen_reales_del_gen(gen_id):
    fen_reales = gene_phenotype_dict[str(gen_id)]
    return set(fen_reales)


def whats_your_set(mph=0.5,iph=0.5,type_of_noise="normal"):
    """
Esta función es para elegir el set simulado de la base de datos, le das un
parámetro de mph e iph y el tipo de distribución de ruido que quieras, y te
devuelve el set
    """
    if type_of_noise == 'normal':
        filename = f'{PATH}data/simulated/{type_of_noise}_simulations/mph_mean_{mph}_mph_std0.1_iph_mean{iph}_iph_std_0.1.txt'
        with open(filename,'r') as file:
            noised_gene_phenotype_dict = json.load(file)

    elif type_of_noise == 'constant':
        filename = f'{PATH}data/simulated/{type_of_noise}_simulations/mph_{mph}_iph_{iph}.json'
        with open(filename,'r') as file:
            noised_gene_phenotype_dict = json.load(file)


    elif type_of_noise == 'random':
        filename = f'{PATH}data/simulated/{type_of_noise}_simulations/random_simulated_data.json'
        with open(filename,'r') as file:
            noised_gene_phenotype_dict = json.load(file)

    elif type_of_noise == 'gold_standard':
        filename = f'{PATH}data/simulated/gene_phenotype_dict.json'
        with open(filename,'r') as file:
            noised_gene_phenotype_dict = json.load(file)

    return noised_gene_phenotype_dict


def fen_observados_con_ruido(gen_id,
        noised_set,n_sample):

    """
mph e iph corresponden a los dos parámetros de ruido: missing phenotypes e
incorrect phenotypes, respectivamente, el que elijamos ahí va a ser tomado de
la base de datos.
n_sample es el número de fenotipos observados que queremos tomar
    """

    fen_observados_total = dict(noised_set)[str(gen_id)]
    try:
        n_fen_observados = np.random.choice(fen_observados_total,n_sample,replace=False)
    except:
        #Si el n_sample es mayor al n total de fen para ese gen, entonces,
        # directamente le ponemos el fen_total de ese gen (que va a ser menor)
        n_fen_observados = fen_observados_total    #replace=False es para que una vez que elijo un fenotipo, no lo vuelva a elegir.

    return set(n_fen_observados)

## }}}



## {{{

db_gene_to_phenotype = f'{PATH}data/simulated/gene_phenotype_dict.json'
with open(db_gene_to_phenotype,'r') as file:
    gene_phenotype_dict = json.load(file)

def lista_de_genes():
    lista = [i for i in gene_phenotype_dict.keys()]
    return lista

def corrida_de_pesos(type_of_noise="normal",mph=0.1,iph=0.1):

    lista = lista_de_genes()
    noised_set = whats_your_set(mph,iph,type_of_noise)
    list_of_df = []
    i=0
    for gen in lista:
        fen_reales = fen_reales_del_gen(gen)
        fen_observados = fen_observados_con_ruido(gen,noised_set,100)

        especificidad = especificidad_del_gen(fen_observados,fen_reales)
        capitalidad = capitalidad_del_gen(fen_observados,fen_reales)
        similaridad = similaridad_del_gen(fen_observados,fen_reales)
        """
        And now we will save the results in a matrix where the columns are
        gen and the rows are especificidad, capitalidad and similaridad
        """
        # Create a DataFrame for each gene and append it to the list
        df = pd.DataFrame([{
            'gen': gen,
            'especificidad': especificidad,
            'capitalidad': capitalidad,
            'similaridad': similaridad
            }])

        list_of_df.append(df)
        print(f'Calculando...{i/4992*100:.2f}%', end='\r')
        i+=1

    # Concatenate all the DataFrames in the list
    result_df = pd.concat(list_of_df, ignore_index=True)

    return result_df

def plot_correlations(df):
    parameters = ['especificidad', 'capitalidad', 'similaridad']
    n = len(parameters)

    # Create a 3x3 grid of subplots
    fig, axs = plt.subplots(n, n, figsize=(15, 15))

    for i in range(n):
        for j in range(n):
            # Scatter plot for each pair of parameters
            axs[i, j].scatter(df[parameters[i]], df[parameters[j]], alpha=0.5)

            # Set the title to the correlation coefficient
            corr_coef = np.corrcoef(df[parameters[i]], df[parameters[j]])[0, 1]
            axs[i, j].set_title(f'Correlation: {corr_coef:.2f}')

            # Set the x and y labels
            axs[i, j].set_xlabel(parameters[i])
            axs[i, j].set_ylabel(parameters[j])

    plt.tight_layout()
    plt.show()

## }}}



# def raw_phen_to_genes(pacient_phen_list,phenotype_to_genes):
    # """
# Esta función toma una lista de fenotipos y devuelve una lista de todos los genes
# candidatos. No pesa los fenotipos, no pesa los genes. Simplemente devuelve la
# unión de los genes asociados a los fenotipos dados.
# Para eso recorre la base de datos phenotype_to_genes.
    # """


# def pacient_gene_score(pacient_genes_list):
    # """
# Esta función toma la lista de raw_phen_to_genes del paciente y le asigna un score
# a los elementos en base a su capitalidad y su especificidad con respecto a
# pacient_phen_list
    # """

# def pacient_phen_score(pacient_phen_list):
    # """
# Esta función toma la lista de fenotipos y devuelve un diccionario (u otra
# clase) con sus scores.
    # """

# def general_gene_score(genes_to_phenotype):
    # """
# Esta función toma el archivo genes_to_phenotype.txt y ...
# acá hay que hacer redes y centralidades
    # """

# def general_phen_score(phenotype_to_genes):
    # """
# Esta función toma el archivo phenoype_to_genes.txt y le asigna un score a todos
# los fenotipos según la relación 1/n° genes que lo causan. De modo que aquél
# fenotipo que tiene una relación 1:1 tiene el máximo score.
    # """



