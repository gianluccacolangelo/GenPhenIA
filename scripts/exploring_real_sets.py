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
dictionary = convert_to_dict(file_path)
# print(dictionary)

def save_to_json(dictionary, json_file_path):
    with open(f'{PATH}data/real/bitgenia.json', 'w') as json_file:
        json.dump(dictionary, json_file, indent=4, sort_keys=True)
save_to_json(dictionary,"bla")
## }}}


## {{{

## }}}
