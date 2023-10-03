## {{{
import set_simulator_utils as simulator
import json
import incorporating_orpha as orpha
import numpy as np
PATH = '/home/brainy/Desktop/1ercuatri2023/Tesis/GenPhenIA/'
import sys
sys.path.insert(0,'/home/brainy/Desktop/1ercuatri2023/Tesis/GenPhenIA/src')
import phen_gen_weight_functions as pgw

## }}}

## {{{
# simulated_data = simulator.genphen_simulator(n_samples_per_gene=1)
# with open(f'{PATH}/data/simulated/simulated_set.json','w') as file:
    # json.dump(simulated_data,file)

# simulator.genphen_simulator(type_of_noise = "normal")


with open(f'{PATH}/data/clinical_cases/bitgenia.json','r') as file:
    simulated_data = json.load(file)

genes = list(simulated_data.keys())

synthetic_phens = {}
for gene in genes:

    synhtetic_phens =simulator.single_disease_simulator(int(gene))
    synthetic_phens[gene] = synhtetic_phens
## }}}


## {{{
result = orpha.v2_model_evaluation(synthetic_phens)

## }}}

## {{{

bitgenia_results = orpha.v2_model_evaluation(simulated_data)
## }}}

##{{{
#TODO now
# [x] tomar distribuciones de pesos empíricos
# [] agregar feature de pesos empíricos a single_disease_simulator
# el algoritmo tiene que primero elegir un peso en base a la distribución, y
# luego elegir de manera random los términos que tengan ese peso para esa enf.
empirical_weights = [i for x in bitgenia_results[2] for i in x]
phen_freq, counts = np.unique(empirical_weights,return_counts=True)

## }}}
## {{{
import matplotlib.pyplot as plt

bitgenia_rankings = [x for x in bitgenia_results[0] if isinstance(x,int)]
rankings = [x for x in result[0] if isinstance(x,int)]

plt.hist(bitgenia_rankings,20)
plt.hist(rankings,20,alpha=.5)

## }}}
