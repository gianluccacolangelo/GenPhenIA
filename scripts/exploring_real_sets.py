"""
En este script se exploraran los datos reales. Por un lado los de bitgenia y
por otro lado los de ClinVar. Se mirará por ejemplo, la cantidad de fenotipos
observados promedio, se mirará también la cantidad de missing phens y la
cantidad de incorrect phens. Y por último, yendo a un nivel más profundo, se
mirará el nivel de vaguedad de un fenotipo. Esto es, a qué distancia
está del fenotipo específico que causa el gen.
"""

## {{{ IMPORTACIONES
import numpy as np
import random
import sys
sys.path.insert(0,'/home/brainy/Desktop/1ercuatri2023/Tesis/GenPhenIA/src')
import linear_model_lab as lml
import phen_gen_weight_functions as pgw
PATH = "/home/brainy/Desktop/1ercuatri2023/Tesis/GenPhenIA/"
import ast
import json
import scienceplots
import matplotlib.pyplot as plt
## }}}

## {{{ Convertir a diccionario el set que nos dieron

def convert_to_dict(file_path):
    with open(file_path, 'r') as f:
        content = f.read()

    # This will give us a list of tuples
    list_of_tuples = ast.literal_eval(content)

    # Now, we will convert it to a dictionary
    # We'll also convert the string representation of list to an actual list using ast.literal_eval
    dictionary = {t[0]: ast.literal_eval(t[1]) for t in list_of_tuples}

    return dictionary

file_path = f"{PATH}data/real/bitgenia_IDs_tov3.sets"
dictionary_bitgenia = convert_to_dict(file_path)
dictionary_clinvar = convert_to_dict(f"{PATH}data/real/clinvar_IDs_tov3.sets")

# print(dictionary)

def save_to_json(dictionary, json_file_path):
    with open(json_file_path, 'w') as json_file:
        json.dump(dictionary, json_file, indent=4, sort_keys=True)
save_to_json(dictionary_bitgenia,f'{PATH}data/real/bitgenia.json')
save_to_json(dictionary_clinvar,f'{PATH}data/real/clinvar.json')
## }}}


## {{{

#TODO empezar a hacer funciones que me digan cuántos mph e iph tienen cada
# uno de los sets

with open(f'{PATH}data/real/bitgenia.json','r') as f:
    bitgenia = json.load(f)

with open(f'{PATH}data/real/clinvar.json','r') as f:
    clinvar = json.load(f)

bitgenia_total_phenotypes = [len(phen_set) for phen_set in bitgenia.values()]
clinvar_total_phenotypes = [len(phen_set) for phen_set in clinvar.values()]

with plt.style.context(['science','ieee','nature']):
    fig, (ax1,ax2) = plt.subplots(1,2)
    ax1.hist(bitgenia_total_phenotypes,bins=20,alpha=.7)
    ax1.set_title("Bitgenia")
    ax1.set_xlabel("Fenotipos observados totales (N)",size=4)
    ax1.set_ylabel("Frecuencia",size=4)
    ax1.axvline(x=np.mean(bitgenia_total_phenotypes),linestyle="--",
            label=f"media $={np.round(np.mean(bitgenia_total_phenotypes))}$")
    ax1.tick_params(axis='both',which='major',labelsize=3)
    ax1.legend(fontsize=5)
    ax1.text(20,60,f"std $= {np.round(np.std(bitgenia_total_phenotypes))}$",fontsize=4)

    ax2.hist(clinvar_total_phenotypes,bins=20,alpha=.7)
    ax2.set_title("Clinvar")
    ax2.set_xlabel("Fenotipos observados totales (N)",size=4)
    ax2.axvline(x=np.mean(clinvar_total_phenotypes),linestyle="--",
            label=f"media $={np.round(np.mean(clinvar_total_phenotypes))}$")
    ax2.tick_params(axis='both',which='major',labelsize=3)
    ax2.legend(fontsize=5)
    ax2.text(25,360,f"std $= {np.round(np.std(clinvar_total_phenotypes))}$",fontsize=4)

## }}}


## {{{


## }}}
