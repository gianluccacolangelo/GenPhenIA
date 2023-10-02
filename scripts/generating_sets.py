## {{{
import set_simulator_utils as simulator
import json
import incorporating_orpha as orpha

PATH = '/home/brainy/Desktop/1ercuatri2023/Tesis/GenPhenIA/'
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

## {{{
import matplotlib.pyplot as plt

bitgenia_rankings = [x for x in bitgenia_results[0] if isinstance(x,int)]
rankings = [x for x in result[0] if isinstance(x,int)]

plt.hist(bitgenia_rankings,20)
plt.hist(rankings,20,alpha=.5)

## }}}
