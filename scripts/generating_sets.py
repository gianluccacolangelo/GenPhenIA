import set_simulator_utils as simulator
import json

missing_phens = [0.1,0.2,0.3,0.4,0.5]
incorrect_phens = [0.1,0.2,0.3,0.4,0.5]
PATH = '/home/brainy/Desktop/1ercuatri2023/Tesis/GenPhenIA/'


for i in missing_phens:
    for j in incorrect_phens:
        simulated_data = simulator.genphen_simulator(i,j)
        with open(f'{PATH}/data/simulated/simulated_set_mph{i}_iph{j}.json','w') as file:
            json.dump(simulated_data, file)

