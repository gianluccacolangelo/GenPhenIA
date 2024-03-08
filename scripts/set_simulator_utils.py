#!/home/brainy/miniconda3/envs/gian/bin/python3
"""
Este script contiene las funciones que me van a permitir generar set simulados
a partir de la base de datos genes_to_phenotype y phenotype_to_genes.
"""

## {{{ Importaciones
import numpy as np
import random
import csv
import json
import scipy.stats as stats
from collections import Counter

PATH = "/home/brainy/Desktop/Tesis/GenPhenIA/"
import sys

sys.path.insert(0, "/home/brainy/Desktop/Tesis/GenPhenIA/src")
import incorporating_orpha as orpha
import phen_gen_weight_functions as pgw

## }}}


## {{{ txt_to_dict

# REEMPLAZAR POR EL PATH DONDE TENES TU REPOSITORIO GenPhenIA
PATH = "/home/brainy/Desktop/Tesis/GenPhenIA/"


def txt_to_dict(filename):
    """
    Esta función toma ya sea genes_to_phenotype.txt o phenotype_to_genes.txt y la
    convierte en un diccionario con la siguiente forma
    {'gene1':['phen1','phen2','phen3'], gene2: ....} almacenándolo en formato .json
    """
    gene_phenotype_dict = {}
    with open(filename, "r") as file:
        reader = csv.reader(file, delimiter="\t")
        next(reader)
        for row in reader:
            gene_id, gene_symbol, hpo_id, hpo_name = row
            if gene_id not in gene_phenotype_dict:
                gene_phenotype_dict[gene_id] = []
            gene_phenotype_dict[gene_id].append(hpo_id)

    with open(f"{PATH}/data/simulated/phenotype_gene_dict.json", "w") as file:
        json.dump(gene_phenotype_dict, file)


## }}}

## {{{ opening dict
with open(f"{PATH}data/simulated/gene_phenotype_dict.json", "r") as file:
    gene_phenotype_dict = json.load(file)
## }}}


## {{{ generate_incorrect_phenotypes
def generate_incorrect_phenotypes(phenotypes, num_incorrect, all_phenotypes):
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

    incorrect_phenotypes = random.sample(phenotypes_not_in_current_gene, num_incorrect)
    return incorrect_phenotypes


## }}}

## {{{ genphen_simulator


def genphen_simulator(
    missing_phens=[0.1, 0.2, 0.3, 0.4, 0.5],
    incorrect_phens=[0.1, 0.2, 0.3, 0.4, 0.5],
    n_samples_per_gene=1,
    type_of_noise="random",
    genphen_db=f"{PATH}/data/simulated/genes_to_phenotype.txt",
):
    """
    Esta función toma la base de datos genes_to_phenotype.txt y genera un set
    simulado con los porcentajes de ruidos de missing_phens e incorrect_phens.
    Para el cual tenemos tres opciones de distribución:

        - constant: Siempre un mismo porcentaje de ruido (e.g. mph=0.2, iph=0.1).
        - normal: Elije entre la lista con una distribución normal.
        - random: Elije entre la lista de manera aleatoria (uniforme).
    """

    all_phenotypes = set()
    for phenlist in gene_phenotype_dict.values():
        for phen in phenlist:
            all_phenotypes.add(phen)

    simulated_data = []
    n = 0

    if type_of_noise == "random":
        for gene, phenotypes in gene_phenotype_dict.items():
            missing_num = int(
                np.round(len(phenotypes) * np.random.choice(missing_phens, 1))
            )
            incorrect_num = int(
                np.round(len(phenotypes) * np.random.choice(incorrect_phens, 1))
            )
            for i in range(n_samples_per_gene):

                # Generate incomplete phenotypes
                missing_phenotypes = random.sample(phenotypes, missing_num)

                # Generate incorrect phenotypes
                incorrect_phenotypes = generate_incorrect_phenotypes(
                    phenotypes, incorrect_num, all_phenotypes
                )

                # Create modified phenotype list
                simulated_phenotypes = [
                    p for p in phenotypes if p not in missing_phenotypes
                ] + incorrect_phenotypes
                random.shuffle(simulated_phenotypes)

                # Add modified phenotype list and associated gene to the simulated dataset
                simulated_data.append((gene, simulated_phenotypes))

            print(f"Generando set: {n/4921*100:.2f}% | Ruido aleatorio", end="\r")
            n += 1
        with open(
            f"{PATH}data/simulated/{type_of_noise}_simulations/random_simulated_data.json",
            "w",
        ) as file:
            json.dump(simulated_data, file)

    elif type_of_noise == "normal":
        # Preguntamos qué media std queremos
        mph_std = float(input("missing_phens desvío estándar: "))
        mph_mean = float(input("missing phens media: "))
        iph_std = float(input("incorrect phens desvío estándar: "))
        iph_mean = float(input("incorrect phens media: "))

        # Aplicamos la media y el desvío estándar a la distribución dada
        mph_dist = stats.norm(mph_mean, mph_std)
        mph_probs = mph_dist.pdf(missing_phens)
        iph_dist = stats.norm(iph_mean, iph_std)
        iph_probs = iph_dist.pdf(incorrect_phens)

        # Normalizamos las probabilidades
        mph_probs /= mph_probs.sum()
        iph_probs /= iph_probs.sum()

        for gene, phenotypes in gene_phenotype_dict.items():

            # Tomamos un porcentaje con la distribución de probabilidad normal
            missing_num = int(
                np.round(len(phenotypes) * np.random.choice(missing_phens, p=mph_probs))
            )
            incorrect_num = int(
                np.round(
                    len(phenotypes) * np.random.choice(incorrect_phens, p=iph_probs)
                )
            )
            for i in range(n_samples_per_gene):

                # Generate incomplete phenotypes
                missing_phenotypes = random.sample(phenotypes, missing_num)

                # Generate incorrect phenotypes
                incorrect_phenotypes = generate_incorrect_phenotypes(
                    phenotypes, incorrect_num, all_phenotypes
                )

                # Create modified phenotype list
                simulated_phenotypes = [
                    p for p in phenotypes if p not in missing_phenotypes
                ] + incorrect_phenotypes
                random.shuffle(simulated_phenotypes)

                # Add modified phenotype list and associated gene to the simulated dataset
                simulated_data.append((gene, simulated_phenotypes))

            print(
                f"Generando set: {n/4921*100:.2f}% | Ruido con distribución normal | Missing phens: media={mph_mean} y std={mph_std} | Incorrect phens: media={iph_mean} y std={iph_std}   ",
                end="\r",
            )
            n += 1
        with open(
            f"{PATH}data/simulated/{type_of_noise}_simulations/mph_mean_{mph_mean}_mph_std{mph_std}_iph_mean{iph_mean}_iph_std_{iph_std}.txt",
            "w",
        ) as file:
            json.dump(simulated_data, file)

    elif type_of_noise == "constant":
        mph_ratio = float(input("missing_phens ratio: "))
        iph_ratio = float(input("incorrect_phens ratio: "))

        for gene, phenotypes in gene_phenotype_dict.items():
            missing_num = int(np.round(len(phenotypes) * mph_ratio))
            incorrect_num = int(np.round(len(phenotypes) * iph_ratio))
            for i in range(n_samples_per_gene):

                # Generate incomplete phenotypes
                missing_phenotypes = random.sample(phenotypes, missing_num)

                # Generate incorrect phenotypes
                incorrect_phenotypes = generate_incorrect_phenotypes(
                    phenotypes, incorrect_num, all_phenotypes
                )

                # Create modified phenotype list
                simulated_phenotypes = [
                    p for p in phenotypes if p not in missing_phenotypes
                ] + incorrect_phenotypes
                random.shuffle(simulated_phenotypes)

                # Add modified phenotype list and associated gene to the simulated dataset
                simulated_data.append((gene, simulated_phenotypes))
            print(
                f"Generando set: {n/4921*100:.2f}%\nRuido constante, mph={mph_ratio},iph={iph_ratio}",
                end="\r",
            )
            n += 1

        with open(
            f"{PATH}data/simulated/{type_of_noise}_simulations/mph_{mph_ratio}_iph_{iph_ratio}.txt",
            "w",
        ) as file:
            json.dump(simulated_data, file)

    return print("Completado")


## }}}

## {{{

with open(f"{PATH}data/clinical_cases/bitgenia.json", "r") as f:
    bitgenia = json.load(f)

with open(f"{PATH}data/clinical_cases/clinvar.json", "r") as f:
    clinvar = json.load(f)

bitgenia_total_phenotypes = [len(phen_set) for phen_set in bitgenia.values()]
clinvar_total_phenotypes = [len(phen_set) for phen_set in clinvar.values()]

values, counts = np.unique(bitgenia_total_phenotypes, return_counts=True)
probabilities = counts / len(bitgenia_total_phenotypes)

dist = dict(zip(values, counts))


def how_many_observed(distribution=dist, probabilities=probabilities):
    """
    Recibe: dist, un diccionario de frecuencias de número de observaciones

    Devuelve: un entero, el número de observaciones que se van a generar

    Este número se genera a partir de una distribución de probabilidad empírica
    (obtenida observando los casos reales) para el número de fenotipos
    observados.

    Otra opción, a chequear en caso de que las cosas no funcionen bien, es
    en lugar de tener una dist de observados netos, tener una distribución de
    observados respecto al total para cada enfermedad.
    """
    values = list(distribution.keys())

    return np.random.choice(values, p=probabilities)


exact_values, exact_counts = (
    np.array(
        [
            0.0,
            0.03571429,
            0.125,
            0.14285714,
            0.2,
            0.23809524,
            0.25,
            0.27272727,
            0.28571429,
            0.3,
            0.33333333,
            0.375,
            0.4,
            0.42857143,
            0.5,
            0.5625,
            0.57142857,
            0.6,
            0.625,
            0.63636364,
            0.64705882,
            0.66666667,
            0.71428571,
            0.73333333,
            0.75,
            0.8,
            0.83333333,
            0.85714286,
            0.875,
            0.88888889,
            0.9,
            0.91666667,
            0.92857143,
            1.0,
        ]
    ),
    np.array(
        [
            31,
            1,
            2,
            1,
            3,
            1,
            6,
            1,
            1,
            1,
            5,
            1,
            4,
            2,
            13,
            1,
            3,
            4,
            4,
            2,
            1,
            8,
            1,
            1,
            7,
            2,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            45,
        ]
    ),
)

exact_probabilities = exact_counts / exact_counts.sum()

inexact_values, inexact_counts = (
    np.array(
        [
            0.25,
            0.625,
            0.72727273,
            0.76470588,
            0.79487179,
            0.8,
            0.83333333,
            0.84848485,
            0.86666667,
            0.9,
            0.90410959,
            0.91549296,
            0.92857143,
            0.93617021,
            0.93902439,
            0.94067797,
            0.94444444,
            0.94642857,
            0.94871795,
            0.94897959,
            0.95,
            0.95081967,
            0.95918367,
            0.95945946,
            0.95977011,
            0.95982143,
            0.96,
            0.96261682,
            0.96551724,
            0.96761134,
            0.97087379,
            0.97169811,
            0.97385621,
            0.97402597,
            0.97457627,
            0.975,
            0.97560976,
            0.97647059,
            0.97777778,
            0.97945205,
            0.97959184,
            0.98039216,
            0.98060942,
            0.98076923,
            0.98095238,
            0.98230088,
            0.98245614,
            0.98369565,
            0.98387097,
            0.98445596,
            0.98545455,
            0.98550725,
            0.9869281,
            0.98695652,
            0.98734177,
            0.9875,
            0.98866856,
            0.98876404,
            0.98901099,
            0.9893617,
            0.98972603,
            0.98976982,
            0.99,
            0.99009901,
            0.99054374,
            0.99090909,
            0.99112426,
            0.99117647,
            0.99122807,
            0.99159664,
            0.99168975,
            0.9919571,
            0.99275362,
            0.99280576,
            0.9929078,
            0.99300699,
            0.99310345,
            0.99342105,
            0.99367089,
            0.99371069,
            0.99382716,
            0.99383984,
            0.99393939,
            0.9939759,
            0.99411765,
            0.99415205,
            0.99425287,
            0.99444444,
            0.99459459,
            0.99473684,
            0.99484536,
            0.99515738,
            0.99516908,
            0.99521531,
            0.99526066,
            0.99593496,
            0.99595142,
            0.996139,
            0.99616858,
            0.99626866,
            0.99656357,
            0.99688474,
            0.99691358,
            0.99696049,
            0.99700599,
            0.99706745,
            0.99733333,
            0.99833611,
            0.9984127,
            1.0,
        ]
    ),
    np.array(
        [
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            3,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            3,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            2,
            1,
            1,
            1,
            45,
        ]
    ),
)

inexact_probabilities = inexact_counts / inexact_counts.sum()

phen_weight_freq, phen_weight_counts = (
    np.array([2.5, 17.0, 54.5, 89.5]),
    np.array([8, 54, 73, 67]),
)
weight_probabilities = phen_weight_counts / phen_weight_counts.sum()
categories = [
    "Very rare (<4-1%)",
    "Occasional (29-5%)",
    "Frequent (79-30%)",
    "Very frequent (99-80%)",
]


def observed_distributions(
    total_observed, exact_prob=exact_probabilities, inexact_prob=inexact_probabilities
):
    obs_exact_proportion = np.random.choice(exact_values, p=exact_prob)
    exact = int(np.round(obs_exact_proportion * total_observed))
    inexactos_netos = total_observed - exact
    inexact_inespecific = int(
        np.round(np.random.choice(inexact_values, p=inexact_prob) * inexactos_netos)
    )
    inexact_error = inexactos_netos - inexact_inespecific

    return (exact, inexact_inespecific, inexact_error)


def weight_distribution(
    total_observed, weight_prob=weight_probabilities, categories=categories
):
    """
    total_observed son los fenotipos exactos totales que se 'observaron'
    weight_prob son las prob de que se observe un fenotipo de cada categoría
    categories son las categorías de peso, que son 4

    DEVUELVE: un dict de categories:counts, que sigen la distribución weight_prob
    """
    obs_weight_proportion = Counter(
        np.random.choice(categories, total_observed, p=weight_prob)
    )
    return obs_weight_proportion


def weight_choosing(phens_dict, weight_proportions):
    """
    phen_dict es el diccionario de fenotipos: peso que obtenemos para cada
    enfermedad
    weight_proportions es la cantidad de c/u de las proporciones de fenotipos
    que debemos sacar de phen_dict

    DEVUELVE: una lista de términos de phen_dict que tienen la distribución de
    pesos weight_proportions
    """
    selected_terms = []
    for category, count in weight_proportions.items():
        terms = [term for term, weight in phens_dict.items() if weight == category]
        if len(terms) >= count:
            selected_terms.extend(np.random.choice(terms, count, replace=False))

    return selected_terms


def simulate_v1(exact_prob=exact_probabilities, inexact_prob=inexact_probabilities):
    total_obs = how_many_observed()
    return observed_distributions(total_obs)


with open(f"{PATH}data/simulated/gene_phenotype_dict.json", "r") as file:
    gold_standard = json.load(file)

with open(f"{PATH}data/simulated/vague_gene_phenotype_dict.json", "r") as f:
    vague_gene_phenotype = json.load(f)

all_phenotypes = np.unique(
    [phen for phens_list in vague_gene_phenotype.values() for phen in phens_list]
)
random.shuffle(all_phenotypes)


def single_disease_simulator(gene_number, with_orpha=True):
    """
    Esta función crea fenotipos observados para un gen, siguiendo las
    distribuciones observadas en los casos reales.
    """
    synthetic_observations = []
    n_exact, n_inexact, n_err = (
        simulate_v1()
    )  # tomamos la dist de fenotipos exact, inexact y err
    exact_pool = set(gold_standard[str(gene_number)])  # lista de términos exact
    vague_pool = set(vague_gene_phenotype[str(gene_number)])
    # lista de términos inexactos:
    inexact_pool = exact_pool.union(vague_pool) - exact_pool.intersection(vague_pool)

    if n_exact > 0:

        gene_symbol_output = orpha.translate_entrez_to_gene(int(gene_number))
        if gene_symbol_output:  # Check if the result is not empty or None
            gene_symbol = gene_symbol_output.split(", ")[0]
        else:
            gene_symbol = " "
            print(f"No gene symbol found for gene number {gene_number}.")
        exact_phens = []
        if with_orpha == True:
            if gene_symbol != "":
                try:
                    disease = pgw.gene_diseases(gene_symbol)[0]
                except IndexError:
                    print(f"No diseases found for gene {gene_symbol}.")
                    disease = None  # or however you want to handle this case

            if disease:  # checking if disease is not None
                orpha_dict_phens = pgw.disease_phens(disease)
                weight_dist = weight_distribution(n_exact)
                exact_phens = weight_choosing(orpha_dict_phens, weight_dist)
        if len(exact_phens) < n_exact:
            exact_pool = exact_pool - set(exact_phens)
            added_exact_phens = np.random.choice(
                list(exact_pool), n_exact - len(exact_phens)
            )
            exact_phens.extend(added_exact_phens)
        synthetic_observations.append(list(exact_phens))

    if n_inexact > 0:
        inexact_phens = np.random.choice(list(inexact_pool), n_inexact)
        synthetic_observations.append(list(inexact_phens))

    if n_err > 0:
        random_gene = np.random.choice(list(gold_standard.keys()))
        err_pool = vague_gene_phenotype[random_gene]
        err_phens = np.random.choice(err_pool, n_err)
        synthetic_observations.append(list(err_phens))

    synthetic_observations = [
        item for sublist in synthetic_observations for item in sublist
    ]

    return synthetic_observations


# de acá sacamos la proporción de exactos e inexactos de los totales
# de los inexactos también sacamos la proporción de vagos y erróneos
# ambos tienen una distribución polarizada en los extremos,
# mostrando que un set puede ser full exactos, o full inexactos con full
# erróneos.
# TODO actualizar las distribuciones y ponerlas en el set simulator.

clinvar_proportions = []
clinvar_vague_proportions = []
for gene in clinvar.keys():
    # tomamos los exactos
    total_exact_phens = gold_standard[gene]
    # tomamos los fenotipos de clinvar
    clinvar_gene_phens = clinvar[gene]
    # nos quedamos con aquellos de clinvar que están en los exactos
    clinvar_exact = [term for term in clinvar_gene_phens if term in total_exact_phens]
    # calculamos la proporción y lo guardamos en la lista
    exact_proportion = len(clinvar_exact) / len(clinvar_gene_phens)
    clinvar_proportions.append(exact_proportion)

    # ahora hacemos lo mismo para los vagos
    total_vague_phens = vague_gene_phenotype[gene]
    # nos quedamos con aquellos de clinvar que no son exactos
    clinvar_inexact = set(clinvar_gene_phens) - set(clinvar_exact)
    # si no es un set vacío
    if len(clinvar_inexact) > 0:
        # nos quedamos con los inexactos que estan en el set de vagos
        clinvar_vague = [term for term in clinvar_inexact if term in total_vague_phens]
        # calculamos la proporción y lo guardamos en la lista
        vague_proportion = len(clinvar_vague) / len(clinvar_inexact)
        clinvar_vague_proportions.append(vague_proportion)


def random_phen_gen(gene):
    """
    Esta función toma x fenotipos al azar de all_phenotypes
    el x debe estar dado por how_many_observed
    """
    total_observerd = how_many_observed()
    phens = np.random.choice(all_phenotypes, total_observerd)
    return list(phens)


# In [37]: np.unique(clinvar_vague_proportions,retur
# ...: n_counts=True)
# Out[37]:
# (array([0.        , 0.06451613, 0.06666667, 0.1111
# 1111, 0.11764706,
# 0.125     , 0.14285714, 0.15      , 0.1666
# 6667, 0.18181818,
# 0.1875    , 0.19047619, 0.2       , 0.2142
# 8571, 0.22222222,
# 0.23076923, 0.23529412, 0.25      , 0.2727
# 2727, 0.28571429,
# 0.3       , 0.3125    , 0.33333333, 0.3636
# 3636, 0.375     ,
# 0.4       , 0.42857143, 0.44444444, 0.5
# , 0.57142857,
# 0.6       , 0.66666667, 0.75      , 0.8
# , 0.83333333,
# 1.        ]),
# array([417,   1,   1,   1,   1,   2,   6,   1,  5
# 1,   1,   2,   1,  18,
# 1,   1,   1,   1,  14,   1,   4,   1,
# 1,  43,   2,   2,  16,
# 3,   2,  65,   1,   5,  15,   3,   4,
# 1,  78]))

# In [38]: ## {{{
# ...: %paste -q

# In [39]: np.unique(clinvar_vague_proportions,retur
# ...: n_counts=True)
# Out[39]:
# (array([0.        , 0.05882353, 0.0625    , 0.0666
# 6667, 0.09090909,
# 0.14285714, 0.15789474, 0.16666667, 0.1818
# 1818, 0.2       ,
# 0.25      , 0.3       , 0.33333333, 0.3636
# 3636, 0.4       ,
# 0.42857143, 0.5       , 0.66666667, 0.75
# , 0.8       ,
# 1.        ]),
# array([61,  1,  1,  1,  1,  4,  1,  2,  1,  4,  7
# ,  1,  8,  1,  2,  1, 10,
# 4,  2,  1, 15]))

# In [40]: np.unique(clinvar_proportions,return_coun
# ...: ts=True)
# Out[40]:
# (array([0.        , 0.06666667, 0.08333333, 0.0909
# 0909, 0.0952381 ,
# 0.11111111, 0.125     , 0.14285714, 0.1666
# 6667, 0.19047619,
# 0.2       , 0.25      , 0.27272727, 0.2857
# 1429, 0.3125    ,
# 0.33333333, 0.36363636, 0.375     , 0.4
# , 0.42857143,
# 0.43333333, 0.45454545, 0.47619048, 0.5
# , 0.55555556,
# 0.57142857, 0.6       , 0.61111111, 0.625
# , 0.64285714,
# 0.66666667, 0.73333333, 0.75      , 0.7777
# 7778, 0.8       ,
# 0.81818182, 0.83333333, 1.        ]),
# array([47,  1,  1,  1,  1,  1,  2,  1,  1,  1,  3
# ,  6,  1,  1,  1,  4,  1,
# 1,  3,  2,  1,  1,  1, 19,  2,  3,  4,  1
# ,  3,  1,  4,  1,  3,  1,
# 2,  1,  1, 30]))

# In [41]: ## {{{
# ...: %paste -q

# In [42]: np.unique(clinvar_proportions,return_coun
# ...: ts=True)
# Out[42]:
# (array([0.        , 0.05555556, 0.11111111, 0.125
# , 0.13043478,
# 0.13888889, 0.14285714, 0.16666667, 0.1851
# 8519, 0.2       ,
# 0.21428571, 0.22222222, 0.25      , 0.2777
# 7778, 0.28571429,
# 0.3       , 0.30769231, 0.3125    , 0.3333
# 3333, 0.35294118,
# 0.35714286, 0.36363636, 0.375     , 0.3846
# 1538, 0.38888889,
# 0.4       , 0.40740741, 0.4137931 , 0.4166
# 6667, 0.42857143,
# 0.44      , 0.44444444, 0.45833333, 0.4615
# 3846, 0.5       ,
# 0.51515152, 0.52631579, 0.52941176, 0.5333
# 3333, 0.53846154,
# 0.54545455, 0.55555556, 0.57142857, 0.5833
# 3333, 0.6       ,
# 0.60869565, 0.61111111, 0.61538462, 0.625
# , 0.63636364,
# 0.64285714, 0.64705882, 0.66666667, 0.6923
# 0769, 0.7       ,
# 0.7037037 , 0.70588235, 0.71428571, 0.7222
# 2222, 0.72727273,
# 0.73913043, 0.75      , 0.77777778, 0.7826
# 087 , 0.8       ,
# 0.82352941, 0.83333333, 0.84615385, 0.8571
# 4286, 0.88888889,
# 0.9       , 0.90909091, 0.92307692, 0.9411
# 7647, 1.        ]),
# array([434,   1,   1,   1,   1,   1,   2,  19,
# 1,   6,   1,   2,   7,
# 1,   6,   1,   3,   2,  37,   1,   1,
# 4,   1,   2,   1,   9,
# 1,   1,   3,   4,   1,   4,   1,   1,  8
# 4,   1,   1,   1,   2,
# 2,   7,   5,   7,   2,  12,   1,   1,
# 1,   4,   3,   1,   1,
# 23,   3,   2,   1,   1,   1,   1,   2,
# 1,  12,   1,   1,   6,
# 1,   6,   1,   2,   1,   1,   1,   1,
# 1, 109]))

## }}}
