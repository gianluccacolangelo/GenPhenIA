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
## }}}


## {{{ txt_to_dict

#REEMPLAZAR POR EL PATH DONDE TENES TU REPOSITORIO GenPhenIA
PATH = '/home/brainy/Desktop/1ercuatri2023/Tesis/GenPhenIA/'


def txt_to_dict(filename):
    """
Esta función toma ya sea genes_to_phenotype.txt o phenotype_to_genes.txt y la
convierte en un diccionario con la siguiente forma
{'gene1':['phen1','phen2','phen3'], gene2: ....} almacenándolo en formato .json
    """
    gene_phenotype_dict = {}
    with open(filename,'r') as file:
        reader = csv.reader(file,delimiter='\t')
        next(reader)
        for row in reader:
            gene_id, gene_symbol, hpo_id, hpo_name = row
            if gene_id not in gene_phenotype_dict:
                gene_phenotype_dict[gene_id] = []
            gene_phenotype_dict[gene_id].append(hpo_id)

    with open(f'{PATH}/data/simulated/gene_phenotype_dict.json','w') as file:
        json.dump(gene_phenotype_dict, file)
## }}}

## {{{ opening dict
with open(f'{PATH}data/simulated/gene_phenotype_dict.json', 'r') as file:
    gene_phenotype_dict = json.load(file)
## }}}

## {{{ generate_incorrect_phenotypes
def generate_incorrect_phenotypes(phenotypes,num_incorrect,all_phenotypes):
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

    incorrect_phenotypes = random.sample(phenotypes_not_in_current_gene,
            num_incorrect)
    return incorrect_phenotypes
## }}}

## {{{ genphen_simulator

def genphen_simulator(missing_phens=[0.1,0.2,0.3,0.4,0.5],
        incorrect_phens=[0.1,0.2,0.3],
        n_samples_per_gene=1,
        type_of_noise = "random",
        genphen_db=f'{PATH}/data/simulated/genes_to_phenotype.txt'):
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
    n=0



    if type_of_noise == "random":
        for gene, phenotypes in gene_phenotype_dict.items():
            missing_num = int(np.round(len(phenotypes) *
                np.random.choice(missing_phens,1)))
            incorrect_num = int(np.round(len(phenotypes) *
                np.random.choice(incorrect_phens,1)))
            for i in range(n_samples_per_gene):

                # Generate incomplete phenotypes
                missing_phenotypes = random.sample(phenotypes, missing_num)

                # Generate incorrect phenotypes
                incorrect_phenotypes = generate_incorrect_phenotypes(phenotypes,
                        incorrect_num,
                        all_phenotypes)

                # Create modified phenotype list
                simulated_phenotypes = [p for p in phenotypes if p not in missing_phenotypes] + incorrect_phenotypes
                random.shuffle(simulated_phenotypes)

                # Add modified phenotype list and associated gene to the simulated dataset
                simulated_data.append((gene, simulated_phenotypes))

            print(f'Generando set: {n/4921*100:.2f}% | Ruido aleatorio',end='\r')
            n+=1
        with open(f'{PATH}data/simulated/{type_of_noise}_simulations/random_simulated_data.json', 'w') as file:
            json.dump(simulated_data, file)


    elif type_of_noise == "normal":
        #Preguntamos qué media std queremos
        mph_std = float(input("missing_phens desvío estándar: "))
        mph_mean = float(input("missing phens media: "))
        iph_std = float(input("incorrect phens desvío estándar: "))
        iph_mean = float(input("incorrect phens media: "))

        #Aplicamos la media y el desvío estándar a la distribución dada
        mph_dist = stats.norm(mph_mean, mph_std)
        mph_probs = mph_dist.pdf(missing_phens)
        iph_dist = stats.norm(iph_mean, iph_std)
        iph_probs = iph_dist.pdf(incorrect_phens)

        #Normalizamos las probabilidades
        mph_probs /= mph_probs.sum()
        iph_probs /= iph_probs.sum()

        for gene, phenotypes in gene_phenotype_dict.items():

            #Tomamos un porcentaje con la distribución de probabilidad normal
            missing_num = int(np.round(len(phenotypes) *
                np.random.choice(missing_phens,p=mph_probs)))
            incorrect_num = int(np.round(len(phenotypes) *
                np.random.choice(incorrect_phens,p=iph_probs)))
            for i in range(n_samples_per_gene):

                # Generate incomplete phenotypes
                missing_phenotypes = random.sample(phenotypes, missing_num)

                # Generate incorrect phenotypes
                incorrect_phenotypes = generate_incorrect_phenotypes(phenotypes,
                        incorrect_num,
                        all_phenotypes)

                # Create modified phenotype list
                simulated_phenotypes = [p for p in phenotypes if p not in missing_phenotypes] + incorrect_phenotypes
                random.shuffle(simulated_phenotypes)

                # Add modified phenotype list and associated gene to the simulated dataset
                simulated_data.append((gene, simulated_phenotypes))

            print(f'Generando set: {n/4921*100:.2f}% | Ruido con distribución normal | Missing phens: media={mph_mean} y std={mph_std} | Incorrect phens: media={iph_mean} y std={iph_std}   ',end="\r")
            n+=1
        with open(f'{PATH}data/simulated/{type_of_noise}_simulations/mph_mean_{mph_mean}_mph_std{mph_std}_iph_mean{iph_mean}_iph_std_{iph_std}.txt','w') as file:
            json.dump(simulated_data,file)


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
                incorrect_phenotypes = generate_incorrect_phenotypes(phenotypes,
                        incorrect_num,
                        all_phenotypes)

                # Create modified phenotype list
                simulated_phenotypes = [p for p in phenotypes if p not in missing_phenotypes] + incorrect_phenotypes
                random.shuffle(simulated_phenotypes)

                # Add modified phenotype list and associated gene to the simulated dataset
                simulated_data.append((gene, simulated_phenotypes))
            print(f'Generando set: {n/4921*100:.2f}%\nRuido constante, mph={mph_ratio},iph={iph_ratio}',
                    end='\r')
            n+=1

        with open(f'{PATH}data/simulated/{type_of_noise}_simulations/mph_{mph_ratio}_iph_{iph_ratio}.txt','w') as file:
            json.dump(simulated_data,file)

    return print("Completado")

## }}}



