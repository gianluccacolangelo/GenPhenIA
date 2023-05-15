## {{{
import set_simulator_utils as simulator
import json

PATH = '/home/brainy/Desktop/1ercuatri2023/Tesis/GenPhenIA/'
## }}}

## {{{
# simulated_data = simulator.genphen_simulator(n_samples_per_gene=1)
# with open(f'{PATH}/data/simulated/simulated_set.json','w') as file:
    # json.dump(simulated_data,file)

simulator.genphen_simulator(type_of_noise = "normal")

## }}}
