## {{{
import set_simulator_utils as simulator
import json
import incorporating_orpha as orpha
import numpy as np

PATH = "/home/brainy/Desktop/Tesis/GenPhenIA/"
import sys

sys.path.insert(0, "/home/brainy/Desktop/Tesis/GenPhenIA/src")
import phen_gen_weight_functions as pgw
import linear_model_lab as lml
import matplotlib.pyplot as plt
import csv

## }}}

## {{{
# simulated_data = simulator.genphen_simulator(n_samples_per_gene=1)
# with open(f'{PATH}/data/simulated/simulated_set.json','w') as file:
# json.dump(simulated_data,file)

# simulator.genphen_simulator(type_of_noise = "normal")


with open(f"{PATH}/data/clinical_cases/bitgenia.json", "r") as file:
    simulated_data = json.load(file)

all_genes = list(simulator.gold_standard.keys())
genes = list(simulated_data.keys())

# synthetic_phens = {}
# for gene in all_genes:
# variation
# for i in range(50):
# synhtetic_phens=simulator.single_disease_simulator(int(gene),with_orpha=False)
# synthetic_phens[gene] = synhtetic_phens

## }}}

## {{{


def get_variation_for_gene(gene, N):
    variations = []
    for _ in range(N):
        variations.append(simulator.single_disease_simulator(gene, with_orpha=False))
    return variations


def get_all_variation_for_all_genes(genes, N):
    data = {}
    count = 0
    total = len(genes)
    for gene in genes:
        print(f"{count}/{total}")
        count += 1
        try:
            data[gene] = get_variation_for_gene(gene, N)
        except:
            continue
    return data


def store_in_csv(data, file_name="simulated_dataset.csv"):
    with open(file_name, "w", newline="") as file:
        writer = csv.writer(file)

        writer.writerow(["Gene", "Phenotypes"])

        # Write data
        for gene, variations in data.items():
            for variation in variations:
                writer.writerow([gene, ", ".join(variation)])


## }}}

## {{{
data = get_all_variation_for_all_genes(all_genes, 100)
store_in_csv(data, file_name="simulated_dataset_100_per_gene_without_orpha.csv")
## }}}


##{{{

synthetic_phens_2 = {}
for gene in genes:

    synhtetic_phens_2 = simulator.single_disease_simulator(int(gene))
    synthetic_phens_2[gene] = synhtetic_phens_2


##}}}


## {{{
# result = orpha.v2_model_evaluation(synthetic_phens)
result_2 = orpha.v2_model_evaluation(synthetic_phens_2)
bitgenia_results = orpha.v2_model_evaluation(simulated_data)
## }}}
## {{{

negative_control_result = orpha.v2_model_evaluation(negative_control_phens)

## }}}

## {{{
import pickle

# with open('simulación_con_orpha_true.pkl', 'wb') as file:
# pickle.dump(result, file)
with open("simulación_con_orpha_false.pkl", "wb") as file:
    pickle.dump(result_2, file)
with open("bitgenia_results.pkl", "wb") as file:
    pickle.dump(bitgenia_results, file)
## }}}

## {{{
import pickle

# with open('simulación_con_orpha_true.pkl', 'rb') as file:
# result = pickle.load(file)
with open("simulación_con_orpha_false.pkl", "rb") as file:
    result_2 = pickle.load(file)
with open("bitgenia_results.pkl", "rb") as file:
    bitgenia_results = pickle.load(file)

## }}}

## {{{

negative_control_rankings = [1000 for _ in range(101)]
bitgenia_rankings = [x for x in bitgenia_results[0] if isinstance(x, int)]
# rankings = [x for x in result[0] if isinstance(x,int)]
rankings_2 = [x for x in result_2[0] if isinstance(x, int)]
# negative_control_rankings = [x for x in negative_control_result[0] if isinstance(x, int)]
## }}}


##{{{


recall_at_N = []

for case in [bitgenia_rankings, rankings_2, negative_control_rankings]:
    model_recall = []
    for i in range(1, 101):
        recall = lml.percent_below_x(case, i)
        model_recall.append(recall)
    recall_at_N.append(model_recall)


with plt.style.context(["science", "ieee", "nature"]):
    fig, ax = plt.subplots()
    ax.plot(range(1, 101), recall_at_N[0], label="casos reales")
    ax.plot(range(1, 101), recall_at_N[2], label="control negativo")
    ax.plot(range(1, 101), recall_at_N[1], label="simulados")
    ax.set_xlabel("TOP\\#")
    ax.set_ylabel("\\% dentro del TOP\\#")
    ax.set_ylim(0, 1)
    ax.set_xlim(1, 100)
    ax.legend(fontsize=5)

##}}}


##{{{
"""
Acá la idea es estudiar la distribución pesos empíricos y distancias entre los
sintéticos y los reales, para ver qué es lo que los está poniendo más arriba
"""

# result_distncias = [y for x in result[1] for y in x]
# result_2_distancias = [y for x in result_2[1] for y in x]
# bitgenia_distancias = [y for x in bitgenia_results[1] for y in x]


# with plt.style.context(['science','ieee','nature']):
# fig,ax = plt.subplots()
# ax.hist(result_distncias,bins=8,label='simulados',alpha=0.5)
# ax.hist(result_2_distancias,bins=8,label='simulados 2',alpha=0.5)
# ax.hist(bitgenia_distancias,bins=5,label='bitgenia',alpha=0.5)
# ax.set_xlabel('Distancia')
# ax.set_ylabel('Frecuencia')
# ax.legend(fontsize=5)


from collections import Counter

result_total_obs = dict(sorted(Counter([y for x in result[1] for y in x]).items()))
result_2_total_obs = dict(sorted(Counter([y for x in result_2[1] for y in x]).items()))
bitgenia_total_obs = dict(
    sorted(Counter([y for x in bitgenia_results[1] for y in x]).items())
)

with plt.style.context(["science", "ieee", "nature"]):
    fig, ax = plt.subplots()
    ax.bar(
        result_total_obs.keys(), result_total_obs.values(), label="simulados", alpha=0.5
    )
    ax.bar(
        result_2_total_obs.keys(),
        result_2_total_obs.values(),
        label="simulados 2",
        alpha=0.5,
    )
    ax.bar(
        bitgenia_total_obs.keys(),
        bitgenia_total_obs.values(),
        label="bitgenia",
        alpha=0.5,
    )
    ax.set_xlabel("Observaciones")
    ax.set_ylabel("Frecuencia")
    ax.legend(fontsize=5)

# with plt.style.context(['science','ieee','nature']):
# fig,ax = plt.subplots()
# ax.hist(result_total_obs,bins=7,label='simulados',alpha=0.5)
# ax.hist(result_2_total_obs,bins=9,label='simulados 2',alpha=0.5)
# ax.hist(bitgenia_total_obs,bins=5,label='bitgenia',alpha=0.5)
# ax.set_xlabel('Observaciones')
# ax.set_ylabel('Frecuencia')
# ax.legend(fontsize=5)

##}}}


## {{{

bitgenia_phens = list(simulated_data.values())
simulated_phens = list(synthetic_phens.values())


def jaccard_sim(x, y):
    return len(x.intersection(y)) / len(x.union(y))


similarities = []
for i in range(0, 159):
    # now we calculate the jaccard similarity between the real and the simulated
    # phenotypes
    # set(bitgenia_phens[i])
    # set(simulated_phens[i])
    sim = jaccard_sim(set(bitgenia_phens[i]), set(simulated_phens[i]))
    similarities.append(sim)

with plt.style.context(["science", "ieee", "nature"]):
    fig, ax = plt.subplots()
    ax.hist(similarities, bins=10)
    ax.set_xlabel("Similaridad de Jaccard", fontsize=6)
    ax.set_ylabel("Frecuencia", fontsize=6)
    ax.set_title(
        "Similaridad entre los \nfenotipos reales y sintéticos\npara el mismo gen",
        fontsize=6,
    )
    ax.set_xlim(0, 1)
## }}}


## {{{

exact_clinv_values, exact_clinv_counts = (
    np.array(
        [
            0.0,
            0.05555556,
            0.11111111,
            0.125,
            0.13043478,
            0.13888889,
            0.14285714,
            0.16666667,
            0.18518519,
            0.2,
            0.21428571,
            0.22222222,
            0.25,
            0.27777778,
            0.28571429,
            0.3,
            0.30769231,
            0.3125,
            0.33333333,
            0.35294118,
            0.35714286,
            0.36363636,
            0.375,
            0.38461538,
            0.38888889,
            0.4,
            0.40740741,
            0.4137931,
            0.41666667,
            0.42857143,
            0.44,
            0.44444444,
            0.45833333,
            0.46153846,
            0.5,
            0.51515152,
            0.52631579,
            0.52941176,
            0.53333333,
            0.53846154,
            0.54545455,
            0.55555556,
            0.57142857,
            0.58333333,
            0.6,
            0.60869565,
            0.61111111,
            0.61538462,
            0.625,
            0.63636364,
            0.64285714,
            0.64705882,
            0.66666667,
            0.69230769,
            0.7,
            0.7037037,
            0.70588235,
            0.71428571,
            0.72222222,
            0.72727273,
            0.73913043,
            0.75,
            0.77777778,
            0.7826087,
            0.8,
            0.82352941,
            0.83333333,
            0.84615385,
            0.85714286,
            0.88888889,
            0.9,
            0.90909091,
            0.92307692,
            0.94117647,
            1.0,
        ]
    ),
    np.array(
        [
            434,
            1,
            1,
            1,
            1,
            1,
            2,
            19,
            1,
            6,
            1,
            2,
            7,
            1,
            6,
            1,
            3,
            2,
            37,
            1,
            1,
            4,
            1,
            2,
            1,
            9,
            1,
            1,
            3,
            4,
            1,
            4,
            1,
            1,
            84,
            1,
            1,
            1,
            2,
            2,
            7,
            5,
            7,
            2,
            12,
            1,
            1,
            1,
            4,
            3,
            1,
            1,
            23,
            3,
            2,
            1,
            1,
            1,
            1,
            2,
            1,
            12,
            1,
            1,
            6,
            1,
            6,
            1,
            2,
            1,
            1,
            1,
            1,
            1,
            109,
        ]
    ),
)


inexact_clinv_values, inexact_clinv_counts = (
    np.array(
        [
            0.0,
            0.05882353,
            0.0625,
            0.06666667,
            0.09090909,
            0.14285714,
            0.15789474,
            0.16666667,
            0.18181818,
            0.2,
            0.25,
            0.3,
            0.33333333,
            0.36363636,
            0.4,
            0.42857143,
            0.5,
            0.66666667,
            0.75,
            0.8,
            1.0,
        ]
    ),
    np.array([61, 1, 1, 1, 1, 4, 1, 2, 1, 4, 7, 1, 8, 1, 2, 1, 10, 4, 2, 1, 15]),
)

bins = np.linspace(inexact_clinv_values.min(), inexact_clinv_values.max(), 15)
hist_counts, _ = np.histogram(
    inexact_clinv_values, bins=bins, weights=inexact_clinv_counts
)

with plt.style.context(["science", "ieee", "nature"]):
    fig, ax = plt.subplots()
    ax.bar(bins[:-1], hist_counts, width=(bins[1] - bins[0]), align="edge", alpha=0.7)
    ax.set_xlabel("Proporción", fontsize=6)
    ax.set_ylabel("Frecuencia")
    ax.set_title(
        "Distribución de los fenotipos\nvagos para los inexactos\nregistrados de ClinVar",
        fontsize=6,
    )
    ax.set_xlim(0, 1)
## }}}
