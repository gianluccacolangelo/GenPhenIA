"""
Este modulo almacena las funciones que le dan peso a fenotipos y genotipos,
tanto en general, como para los pacientes en particular.
"""

## {{{ IMPORTACIONES
import numpy as np
import json
import pandas as pd
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
from lxml import etree
import sys

sys.path.insert(0, "/home/brainy/Desktop/Tesis/GenPhenIA/scripts/")
import incorporating_orpha as orpha

PATH = "/home/brainy/Desktop/Tesis/GenPhenIA/"
## }}}


## {{{ funciones de peso de genes candidatos

with open(f"{PATH}config/phen_properties.csv", "r") as file:
    phen_promiscuity_dict = pd.read_csv(file, comment="#", sep="\t")
    # weight_dict = phen_promiscuity_dict.set_index('Unnamed: 0')['c/(a+c)'].to_dict()
    weight_dict = phen_promiscuity_dict.set_index("Unnamed: 0")[
        "frecuencia relativa"
    ].to_dict()

# link tree abre el xml que asocia cada orphacode con un mimcode
link_tree = ET.parse(f"{PATH}data/ORPHA/alignments_omim_orpha.xml")
link_root = link_tree.getroot()
mim_orpha_dict = {}

for disorder in link_root.findall(".//Disorder"):

    ext_ref_list = disorder.find("ExternalReferenceList")
    orphacode = disorder.find("OrphaCode").text

    if ext_ref_list is not None:
        # Iterate over all 'ExternalReference' elements in the list
        for ext_ref in ext_ref_list.findall("ExternalReference"):
            # Extract and print information from each 'ExternalReference'
            source = (
                ext_ref.find("Source").text
                if ext_ref.find("Source") is not None
                else None
            )
            reference = (
                ext_ref.find("Reference").text
                if ext_ref.find("Reference") is not None
                else None
            )

            if source == "OMIM":
                mim_orpha_dict[reference] = [orphacode]

gene_disease_tree = ET.parse(f"{PATH}data/ORPHA/gene_disease.xml")
gene_disease_root = gene_disease_tree.getroot()

phen_disease_tree = ET.parse(f"{PATH}data/ORPHA/phen_disease.xml")
phen_disease_root = phen_disease_tree.getroot()

epidemiology_tree = ET.parse(f"{PATH}data/ORPHA/epidemiology.xml")
epidemiology_root = epidemiology_tree.getroot()

classification_tree = ET.parse(
    f"{PATH}data/ORPHA/rare_genetic_diseases_classification.xml"
)
classification_root = classification_tree.getroot()


natural_history_tree = ET.parse(
    f"{PATH}data/ORPHA/natural_history_of_rare_diseases.xml"
)
natural_history_root = natural_history_tree.getroot()


def especificidad_del_gen(fenotipos_observados, fenotipos_del_gen):
    """
    fenotipos_observados y fenotipos_del_gen tienen que ser sets (conjuntos), esta
    función devolverá la fracción de fenotiops_del_gen que están en fenotipos_observados

    Recordar que siempre hablamos del gen candidato para el conj de fenotipos obs.
    """
    j = fenotipos_observados.intersection(fenotipos_del_gen)
    i = set(fenotipos_observados) - j
    k = set(fenotipos_del_gen) - j
    # if type_of_func == "log":
    # result = len(j)/(np.log10(1+len(k))+len(j))
    # elif type_of_func == "exp":
    # result = len(j)/(10**(len(k))+len(j))
    # elif type_of_func == "lineal":
    # result = len(j)/(len(k)+len(j))
    return len(j) / (np.log10(len(k) + 1) + len(j))


def capitalidad_del_gen(fenotipos_observados, fenotipos_del_gen):
    """
    esta función fevolverá la fracción de fenotiops_observados que están en fenotipos_del_gen

    Recordar que siempre hablamos del gen candidato para el conj de fenotipos obs.
    """
    j = fenotipos_observados.intersection(fenotipos_del_gen)
    i = len(set(fenotipos_observados) - j)
    k = len(set(fenotipos_del_gen) - j)
    j = len(j)
    return j / (i + j)


def similaridad_del_gen(fenotipos_observados, fenotipos_del_gen):
    """
    Esta función devolverá la intersección de fenotipos observados y fenotipos del
    gen sobre la unión de ambos.

    Recordar que siempre hablamos del gen candidato para el conj de fenotipos obs.
    """
    j = fenotipos_observados.intersection(fenotipos_del_gen)
    i = len(set(fenotipos_observados) - j)
    k = len(set(fenotipos_del_gen) - j)
    j = len(j)
    return j / (j + i + np.log10(1 + k))


def parametro(
    fenotipos_observados, fenotipos_del_gen, nuevo_parametro, alpha=1, beta=1, gamma=1
):
    """
    Esta función está para probar nuevas métricas, donde nuevo_parametro es un
    entero de la nueva métrica
    """

    j = fenotipos_observados.intersection(fenotipos_del_gen)
    i = set(fenotipos_observados) - j
    k = set(fenotipos_del_gen) - j

    weighted_j = sum([weight_dict.get(phen, 0) for phen in j])
    weighted_i = sum([weight_dict.get(phen, 0) for phen in i])
    weighted_k = sum([weight_dict.get(phen, 0) for phen in k])

    if nuevo_parametro == 0:
        result = alpha * weighted_j - beta * weighted_i - gamma * weighted_k
    elif nuevo_parametro == 1:
        result = alpha * len(j) - beta * len(i) - gamma * len(k)
    elif nuevo_parametro == 2:
        result = weighted_j - weighted_i
    elif nuevo_parametro == 3:
        result = weighted_j / (1 + weighted_i + weighted_k)
    elif nuevo_parametro == 4:
        result = weighted_j
    elif nuevo_parametro == 5:
        result = -weighted_i
    elif nuevo_parametro == 6:
        result = -weighted_k
    elif nuevo_parametro == 7:
        result = -weighted_i - weighted_k
    return result


## }}}

## {{{ funciones de peso de fenotipos


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
db = f"{PATH}data/simulated/"
with open(f"{db}gene_phenotype_dict.json", "r") as file:
    gene_phenotype_dict = json.load(file)
with open(f"{db}vague_gene_phenotype_dict.json", "r") as file:
    vague_gene_phenotype_dict = json.load(file)


def fen_reales_del_gen(gen_id, precision):
    if precision == True:
        fen_reales = gene_phenotype_dict[str(gen_id)]
    elif precision == False:
        fen_reales = vague_gene_phenotype_dict[str(gen_id)]
    return set(fen_reales)


def whats_your_set(mph=0.5, iph=0.5, type_of_noise="normal"):
    """
    Esta función es para elegir el set simulado de la base de datos, le das un
    parámetro de mph e iph y el tipo de distribución de ruido que quieras, y te
    devuelve el set
    """
    if type_of_noise == "normal":
        filename = f"{PATH}data/simulated/{type_of_noise}_simulations/mph_mean_{mph}_mph_std0.1_iph_mean{iph}_iph_std_0.1.txt"
        with open(filename, "r") as file:
            noised_gene_phenotype_dict = json.load(file)

    elif type_of_noise == "constant":
        filename = (
            f"{PATH}data/simulated/{type_of_noise}_simulations/mph_{mph}_iph_{iph}.json"
        )
        with open(filename, "r") as file:
            noised_gene_phenotype_dict = json.load(file)

    elif type_of_noise == "random":
        filename = f"{PATH}data/simulated/{type_of_noise}_simulations/random_simulated_data.json"
        with open(filename, "r") as file:
            noised_gene_phenotype_dict = json.load(file)

    elif type_of_noise == "gold_standard":
        filename = f"{PATH}data/simulated/gene_phenotype_dict.json"
        with open(filename, "r") as file:
            noised_gene_phenotype_dict = json.load(file)

    return noised_gene_phenotype_dict


def fen_observados_con_ruido(gen_id, noised_set, n_sample):
    """
    mph e iph corresponden a los dos parámetros de ruido: missing phenotypes e
    incorrect phenotypes, respectivamente, el que elijamos ahí va a ser tomado de
    la base de datos.
    n_sample es el número de fenotipos observados que queremos tomar

    DEVUELVE: un set de fenotipos

    """

    fen_observados_total = dict(noised_set)[str(gen_id)]
    try:
        n_fen_observados = np.random.choice(
            fen_observados_total, n_sample, replace=False
        )
    except:
        # Si el n_sample es mayor al n total de fen para ese gen, entonces,
        # directamente le ponemos el fen_total de ese gen (que va a ser menor)
        n_fen_observados = fen_observados_total  # replace=False es para que una vez que elijo un fenotipo, no lo vuelva a elegir.

    return set(n_fen_observados)


parent_map = {c: p for p in gene_disease_root.iter() for c in p}


def gene_diseases(gene_symbol):
    """
    A esta función le damos un gen y nos devuelve la lista de enfermedades
    asociadas
    """
    # Find all 'Gene' elements with the given symbol
    genes = gene_disease_root.findall(f".//Gene[Symbol='{gene_symbol}']")

    # Initialize a list to store the associated diseases
    diseases = []

    # Iterate over all found 'Gene' elements
    for gene in genes:
        # Navigate up the tree structure to find the 'Name' of the disease
        disorder_gene_association = parent_map[gene]
        disorder_gene_association_list = parent_map[disorder_gene_association]
        disorder = parent_map[disorder_gene_association_list]
        disease_name = disorder.find("OrphaCode").text

        # Add the disease name to the list
        diseases.append(disease_name)

    return diseases


def disease_genes(disease):
    """
    A esta función le damos una enfermedad y nos devuelve la lista de genes
    asociados
    """

    disorders = gene_disease_root.findall(f".//Disorder[OrphaCode='{disease}']")

    # Initialize a list to store the associated genes
    genes = []

    # Iterate over all found 'Disorder' elements
    for disorder in disorders:
        # Find all 'Gene' elements associated with the disorder
        disorder_gene_associations = disorder.findall(".//DisorderGeneAssociation")
        for disorder_gene_association in disorder_gene_associations:
            gene = disorder_gene_association.find(".//Gene")
            if gene is not None:
                # Get the 'Symbol' of the gene
                gene_symbol = gene.find("Symbol").text
                # Add the gene symbol to the list
                genes.append(gene_symbol)

    return genes


def disease_phens(disease, weight=True):
    """
    A esta función le damos una enfermedad y nos devuelve la lista de fenotipos
    pesados.
    """

    # Find all 'Disorder' elements with the given disease name
    disorders = phen_disease_root.findall(f".//Disorder[OrphaCode='{disease}']")

    # Initialize a list to store the associated phenotypes
    phenotypes = []

    # Iterate over all found 'Disorder' elements
    for disorder in disorders:
        # Find all 'HPODisorderAssociation' elements associated with the disorder
        hpo_disorder_associations = disorder.findall(".//HPODisorderAssociation")
        for hpo_disorder_association in hpo_disorder_associations:
            # Get the 'HPOTerm' and 'HPOFrequency' of the phenotype
            hpo_term = hpo_disorder_association.find(".//HPOId").text
            hpo_frequency = (
                hpo_disorder_association.find(".//HPOFrequency/Name").text
                if weight
                else None
            )
            hpo_diagnostic_criteria = hpo_disorder_association.find(
                ".//DiagnosticCriteria"
            ).text
            # Add the phenotype and its frequency to the list
            phenotypes.append((hpo_term, hpo_frequency))

    return dict(phenotypes)


class_parent_map = {c: p for p in classification_root.iter() for c in p}


# obtenemos los 32 elementos de clasificaciones en una lista
l1_classifications_elements = classification_root.find(
    ".//ClassificationNodeChildList"
).findall("ClassificationNode")

l2_classifications_elements = []
l3_classifications_elements = []

for node in l1_classifications_elements:
    # Level 2: Find all 'ClassificationNode' elements within the current 'ClassificationNode'
    l2_classifications_elements.extend(
        node.findall(".//ClassificationNodeChildList/ClassificationNode")
    )

for node in l2_classifications_elements:
    # Level 3: Find all 'ClassificationNode' elements within the current 'ClassificationNode'
    l3_classifications_elements.extend(
        node.findall(".//ClassificationNodeChildList/ClassificationNode")
    )


def disease_classification(orphacode, level=1):
    """
    A esta función le damos un orphacode y nos devuelve su clasificación en el
    nivel que le pidamos. El 0 corresponde a 'Rare genetic disease' y el segundo a
    32 clasificaciones diferentes.
    """

    disorder = classification_root.find(f".//Disorder[OrphaCode='{orphacode}']")

    if disorder is None:
        return None

    node = disorder

    # usamos un while para que siga iterando hasta que llegue al nivel que nos
    # interesa
    if level == 1:
        classifications_elements = l1_classifications_elements
    elif level == 2:
        classifications_elements = l2_classifications_elements
    elif level == 3:
        classifications_elements = l3_classifications_elements

    i = 0
    while node not in classifications_elements:
        if i == 20:
            return None
        node = class_parent_map.get(node)
        i += 1

    return (node.find(".//Name").text, node.find(".//OrphaCode").text)


def get_average_age_of_onset(gene_code):
    gene_symbol = orpha.translate_entrez_to_gene(gene_code).split(", ")[0]
    orpha_code = int(gene_diseases(gene_symbol)[0])

    # Iterate over all disorders
    for disorder in natural_history_root.find("DisorderList"):
        # Check if the orpha code matches
        if disorder.find("OrphaCode").text == str(orpha_code):
            # Get the list of average age of onset
            average_age_of_onset_list = disorder.find("AverageAgeOfOnsetList")
            if average_age_of_onset_list is not None:
                # Extract all the average age of onset
                ages_of_onset = [
                    average_age_of_onset.find("Name").text.strip()
                    for average_age_of_onset in average_age_of_onset_list
                ]
                # Filter out empty strings
                ages_of_onset = list(filter(None, ages_of_onset))
                # Return the list of average age of onset
                return ages_of_onset
    # If no disorder with the specified orpha code was found, return None
    return None


def phen_diseases(phen_id):
    """
    Esta función devuelve todas las enfermedades asociadas a al menos 1 fenotipo
    de una lista
    """


## }}}


## {{{

db_gene_to_phenotype = f"{PATH}data/simulated/gene_phenotype_dict.json"
with open(db_gene_to_phenotype, "r") as file:
    gene_phenotype_dict = json.load(file)


def lista_de_genes():
    lista = [i for i in gene_phenotype_dict.keys()]
    return lista


def corrida_de_pesos(type_of_noise="normal", mph=0.1, iph=0.1):

    lista = lista_de_genes()
    noised_set = whats_your_set(mph, iph, type_of_noise)
    list_of_df = []
    i = 0
    for gen in lista:
        fen_reales = fen_reales_del_gen(gen)
        fen_observados = fen_observados_con_ruido(gen, noised_set, 100)

        especificidad = especificidad_del_gen(fen_observados, fen_reales)
        capitalidad = capitalidad_del_gen(fen_observados, fen_reales)
        similaridad = similaridad_del_gen(fen_observados, fen_reales)
        """
        And now we will save the results in a matrix where the columns are
        gen and the rows are especificidad, capitalidad and similaridad
        """
        # Create a DataFrame for each gene and append it to the list
        df = pd.DataFrame(
            [
                {
                    "gen": gen,
                    "especificidad": especificidad,
                    "capitalidad": capitalidad,
                    "similaridad": similaridad,
                }
            ]
        )

        list_of_df.append(df)
        print(f"Calculando...{i/4992*100:.2f}%", end="\r")
        i += 1

    # Concatenate all the DataFrames in the list
    result_df = pd.concat(list_of_df, ignore_index=True)

    return result_df


def plot_correlations(df):
    parameters = ["especificidad", "capitalidad", "similaridad"]
    n = len(parameters)

    # Create a 3x3 grid of subplots
    fig, axs = plt.subplots(n, n, figsize=(15, 15))

    for i in range(n):
        for j in range(n):
            # Scatter plot for each pair of parameters
            axs[i, j].scatter(df[parameters[i]], df[parameters[j]], alpha=0.5)

            # Set the title to the correlation coefficient
            corr_coef = np.corrcoef(df[parameters[i]], df[parameters[j]])[0, 1]
            axs[i, j].set_title(f"Correlation: {corr_coef:.2f}")

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
