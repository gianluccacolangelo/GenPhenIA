import set_simulator_utils as simulator
import json
import numpy as np
import random
from scipy.stats import norm

missing_phens = np.array([0.1,0.2,0.3,0.4,0.5])
incorrect_phens = np.array([0.1,0.2,0.3])
PATH = '/home/brainy/Desktop/1ercuatri2023/Tesis/GenPhenIA/'

def generate_simulated_set(mph_values,
        iph_values,
        mph_mean = 0.3,
        mph_std = 0.1,
        iph_mean = 0.2,
        iph_std = 0.1):
    """
Esta función toma una lista de valores posibles de missing phenotypes e incorrect
phenotypes, con una media y un desvío estándar para cada uno, y genera un set
simulado con una distribución normal de esos errores.
    """

    # Calculate probabilities for mph
    mph_probabilities = norm.pdf(mph_values, mph_mean, mph_std)
    mph_probabilities /= mph_probabilities.sum()  # Normalize so probabilities sum to 1
    # Calculate probabilities for iph
    iph_probabilities = norm.pdf(iph_values, iph_mean, iph_std)
    iph_probabilities /= iph_probabilities.sum()  # Normalize so probabilities sum to 1

    mph = np.random.choice(mph_values,p=mph_probabilities)
    iph = np.random.choice(iph_values,p=iph_probabilities)

    simulated_data = simulator.genphen_simulator(mph,
            iph,
            n_samples_per_gene=1)
    with open(f'{PATH}/data/simulated/simulated_set_mph{mph_mean}_iph{iph_mean}.json','w') as file:
        kkkkkkkkkkkkkkkkkkkkkkkkkkkkk
