"""
En este archivo están los plots de todos los gráficos.
"""

## {{{ importaciones
import sys
sys.path.insert(0,'/home/brainy/Desktop/1ercuatri2023/Tesis/GenPhenIA/src')
sys.path.insert(0,'/home/brainy/Desktop/1ercuatri2023/Tesis/GenPhenIA/scripts')
import incorporating_orpha as orpha
import phen_gen_weight_functions as pgw
import linear_model_lab as lml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import json
import random
PATH = "/home/brainy/Desktop/1ercuatri2023/Tesis/GenPhenIA/"
## }}}



## {{{ CORRIENDO MODELO PARA alpha=beta=gamma=0.33
"""
Esto corre los datos para cada tipo de ruido y calcula el accuracy para
cada N
"""

results = {}
for type_of_noise in ['normal','constant','random','gold_standard']:
    if type_of_noise == 'gold_standard':
        mph_iph_metrics = []
        for fen_samples in range(15):
            list_of_gen_rankings = lml.model_evaluating(0.1,0.1,type_of_noise,fen_samples+1,100)
            mph_iph_metrics.append(lml.percent_below_x(list_of_gen_rankings,10))
        results[f"clean_set"] = mph_iph_metrics


    elif type_of_noise == 'normal' or type_of_noise == 'constant':
        for mph in [0.1,0.5]:
            for iph in [0.1,0.5]:
                mph_iph_metrics = []
                for fen_samples in range(15):
                    list_of_gen_rankings = lml.model_evaluating(mph,iph,type_of_noise,fen_samples+1,100)
                    mph_iph_metrics.append(lml.percent_below_x(list_of_gen_rankings,10))
                results[f"{type_of_noise}: mph={mph}, iph={iph}"] = mph_iph_metrics

    elif type_of_noise == 'random':
        mph_iph_metrics = []
        for fen_samples in range(15):
            list_of_gen_rankings = lml.model_evaluating(0.1,0.1,type_of_noise,fen_samples+1,100)
            mph_iph_metrics.append(lml.percent_below_x(list_of_gen_rankings,10))
        results[f"random"] = mph_iph_metrics


top_10_metrics = pd.DataFrame(results)

## {{{ PLOT
"""
Y acá corremos el plot de top_10_metrics
"""

with open(f"{PATH}output/top_10_metrics.csv", "r") as f:
     top_10_metrics = pd.read_csv(f)

x = np.linspace(0, 10, 100)

# Create a set of line styles to use (you can use any suitable line styles)
styles = ['-', '--', '-.', ':']

# Create a colormap
cmap = plt.get_cmap('tab10')

with plt.style.context(['science','ieee','nature']):
    fig, (ax1,ax2) = plt.subplots(1,2)
    plt.rcParams["text.usetex"] = True
    for i, label in enumerate(top_10_metrics.columns[1:5]):
        ax1.plot(np.arange(1,16), top_10_metrics[label], label=label, color=cmap(i), linestyle=styles[i % len(styles)])
    ax1.plot(np.arange(1,16), top_10_metrics['clean_set'], label='clean_set',
            color='black')
    ax1.legend(loc='lower right', fontsize=4)
    ax1.set_xlabel('Total observed phenotypes (N)', fontsize=4)
    ax1.set_ylabel('Accuracy',fontsize=4)
    ax1.set_title('Accuracy of the model for different noise levels', fontsize=4)
    ax1.tick_params(axis='both', which='major', labelsize=3)
    ax1.axhline(y=0.9, color='green', linestyle='--', linewidth=1,alpha=0.5)
    ax1.axvline(x=5, color='green', linestyle='--', linewidth=1,alpha=0.5)
    ax1.set_ylim(0,1)
    plt.subplots_adjust(top=0.9,bottom=0.15,
            left=0.176,right=0.952,
            hspace=0.,wspace=0.1)
    for i, label in enumerate(top_10_metrics.columns[5:9]):
        ax2.plot(np.arange(1,16), top_10_metrics[label], label=label, color=cmap(i), linestyle=styles[i % len(styles)])
    ax2.plot(np.arange(1,16), top_10_metrics['clean_set'], label='clean_set',
            color='black')
    ax2.legend(loc='lower right', fontsize=4)
    ax2.set_xlabel('Total observed phenotypes (N)', fontsize=4)
    # ax2.set_ylabel('Accuracy',fontsize=4)
    ax2.set_title('Accuracy of the model for different noise levels', fontsize=4)
    ax2.tick_params(axis='both', which='major', labelsize=3)
    ax2.axhline(y=0.9, color='green', linestyle='--', linewidth=1,alpha=0.5)
    ax2.axvline(x=5, color='green', linestyle='--', linewidth=1,alpha=0.5)
    ax2.set_ylim(0,1)

plt.show()
## }}}

## }}}



## {{{ CORRIENDO MODELO PARA cada métrica separada

for alpha,beta,gamma in [(1,0,0),(0,1,0),(0,0,1)]:
    results = {}
#alpha -> especificidad
#beta -> capitalidad
#gamma -> similaridad
    for type_of_noise in ['constant','gold_standard']:
        if type_of_noise == 'gold_standard':
            mph_iph_metrics = []
            for fen_samples in range(10):
                list_of_gen_rankings = lml.model_evaluating(0.1,0.1,type_of_noise,
                            fen_samples+1,
                            500,
                            alpha,
                            beta,
                            gamma)
                mph_iph_metrics.append(lml.percent_below_x(list_of_gen_rankings,10))
            results[f"clean_set"] = mph_iph_metrics


        elif type_of_noise == 'normal' or type_of_noise == 'constant':
            for mph in [0.1]:
                for iph in [0.1]:
                    mph_iph_metrics = []
                    for fen_samples in range(10):
                        list_of_gen_rankings = lml.model_evaluating(mph,iph,
                                    type_of_noise,
                                    fen_samples+1,
                                    500,
                                    alpha,
                                    beta,
                                    gamma)
                        mph_iph_metrics.append(lml.percent_below_x(list_of_gen_rankings,10))
                    results[f"{type_of_noise}: mph={mph}, iph={iph}"] = mph_iph_metrics

        elif type_of_noise == 'random':
            mph_iph_metrics = []
            for fen_samples in range(15):
                list_of_gen_rankings = lml.model_evaluating(0.1,0.1,type_of_noise,fen_samples+1,500)
                mph_iph_metrics.append(lml.percent_below_x(list_of_gen_rankings,10))
            results[f"random"] = mph_iph_metrics


        top_10_metrics = pd.DataFrame(results)

        with open(f"{PATH}output/top_10_metrics_{alpha}_{beta}_{gamma}.csv", "w") as f:
            top_10_metrics.to_csv(f)

## {{{ PLOT
with open(f"{PATH}output/top_10_metrics_1_1_0_no_1.csv", "r") as f:
    top_10_metrics = pd.read_csv(f)

x = np.linspace(0, 10, 100)

# Create a set of line styles to use (you can use any suitable line styles)
styles = ['-', '--', '-.', ':']

# Create a colormap
cmap = plt.get_cmap('tab10')

with plt.style.context(['science','ieee','nature']):
    fig, ax1 = plt.subplots()
    plt.rcParams["text.usetex"] = True
    for i, label in enumerate(top_10_metrics.columns[1:2]):
        ax1.plot(np.arange(1,11), top_10_metrics[label], label=label, color=cmap(i), linestyle=styles[i % len(styles)])
    ax1.plot(np.arange(1,11), top_10_metrics['clean_set'], label='clean_set',
            color='black')
    ax1.legend(loc='lower right', fontsize=5)
    ax1.text(5, 0.8, 'acum. accuracy$=0.1999$',fontsize=4)
    ax1.set_xlabel('Total observed phenotypes (N)', fontsize=7)
    ax1.set_ylabel('Accuracy',fontsize=7)
    ax1.set_title('$\\frac{j}{j+10^k}$', fontsize=7)
    ax1.tick_params(axis='both', which='major', labelsize=5)
    # ax1.axhline(y=0.9, color='green', linestyle='--', linewidth=1,alpha=0.5)
    # ax1.axvline(x=5, color='green', linestyle='--', linewidth=1,alpha=0.5)
    ax1.set_ylim(0,1)
    plt.subplots_adjust(top=0.9,bottom=0.15,
            left=0.176,right=0.952,
            hspace=0.,wspace=0.1)
## }}}

## }}}



## {{{
def accuracy(type_of_noise,alpha,beta,gamma,nueva_metrica,n_metrica):
    results = {}
    for type_of_noise in ['constant','gold_standard']:
        if type_of_noise == 'gold_standard':
            mph_iph_metrics = []
            for fen_samples in range(10):
                list_of_gen_rankings = lml.model_evaluating(0.1,0.1,type_of_noise,
                            fen_samples+1,
                            500,
                            alpha,
                            beta,
                            gamma,
                            nueva_metrica=nueva_metrica,
                            n_metrica=n_metrica)
                mph_iph_metrics.append(lml.percent_below_x(list_of_gen_rankings,10))
            results[f"clean_set"] = mph_iph_metrics

        elif type_of_noise == 'constant':
            for mph in [0.1]:
                for iph in [0.1]:
                    mph_iph_metrics = []
                    for fen_samples in range(10):
                        list_of_gen_rankings = lml.model_evaluating(mph,iph,
                                    type_of_noise,
                                    fen_samples+1,
                                    500,
                                    alpha,
                                    beta,
                                    gamma,
                                    nueva_metrica=nueva_metrica,
                                    n_metrica=n_metrica)
                        mph_iph_metrics.append(lml.percent_below_x(list_of_gen_rankings,10))
            results[f"{type_of_noise}: mph={mph}, iph={iph}"] = mph_iph_metrics

    top_10_metrics = pd.DataFrame(results)

    with open(f"{PATH}output/top_10_metrics_{alpha}_{beta}_{gamma}_{nueva_metrica}_{n_metrica}.csv", "w") as f:
        top_10_metrics.to_csv(f)

## }}}


## {{{ Calculando el accuracy acumulado total

"""
Acá vamos a calcular el accuracy acumulado total, para elegir las tres métricas
que mejor se desempeñen para luego probar de a pares.
"""
with open(f'{PATH}output/top_10_metrics_0_1_1_no_3.csv', 'r') as f:
    top_10_metrics = pd.read_csv(f)


def acumulated_accuracy(top_10_metrics):
    """
Esta función calcula el accuracy acumulado a lo largo de N de cada métrica y
devuelve el total promediadio por cada una de las columnas (gold standard y con
error)
    """
    # Calculate sum of each column
    column_sums = top_10_metrics.sum()[1:] #para que no tome la primera columna

    # Calculate mean of sums
    mean_of_sums = column_sums.mean()
    print(column_sums)

    return mean_of_sums/10 #sobre el total de N corridos, para normalizarlo
acumulated_accuracy(top_10_metrics)

:q

## }}}

## de las métricas nuevas, el accuracy acumulado es:

# 1 = 0.3143  j-i-k
# 2 = 0.8362  j-i
# 3 = 0.6779  j/(1+i+k)
# 4 = 0.8319  j #Queda
# 5 = 0.8251  -i
# 6 = 0.1694  -k
# 7 = 0.2442  -i-k

# sim_lin = 0.6663
# sim_log = 0.9054 #Queda
# sim_exp = 0.2023

# cap_lin = 0.8333 #Queda
# cap_log = 0.8321
# cap_exp = 0.8272

# esp_lin = 0.4791
# esp_log = 0.8606 #Queda
# esp_exp = 0.1999


# Quizás capitalidad podemos eliminar y quedarnos solo con j. Luego similaridad
# y especificidad agregan otro tipo de información -i es lo mismo que j

# esp + cap = 0.9018
# esp + sim = 0.9022
# cap + sim = 0.9060

# esp + cap + sim = 0.9046
## Después quedaría calcular también las de similaridad capitalidad y
## especificidad





## {{{ hardcoding
db = f'{PATH}data/real/clinvar.json'

with open(db,'r') as f:
    db_real_clinvar = json.load(f)

db = f'{PATH}data/real/bitgenia.json'

with open(db,'r') as f:
    db_real_bitgenia = json.load(f)

bitgenia_metrics = lml.model_real_evaluating(db_real_bitgenia,1,1,1,"no",1,precision=False)
esp_bitgenia_metrics = lml.model_real_evaluating(db_real_bitgenia,1,0,0,"no",1,precision=False)
cap_bitgenia_metrics = lml.model_real_evaluating(db_real_bitgenia,0,1,0,"no",1,precision=False)
sim_bitgenia_metrics = lml.model_real_evaluating(db_real_bitgenia,0,0,1,"no",1,precision=False)
j_bitgenia_metrics = lml.model_real_evaluating(db_real_bitgenia,0,0,0,"si",4,precision=False)


# clinvar_metrics = lml.model_real_evaluating(db_real_clinvar,1,1,1,"no",1,precision=False)

bitgenia_metrics_setpreciso = lml.model_real_evaluating(db_real_bitgenia,1,1,1,"no",1,precision=True)
esp_bitgenia_metrics_setpreciso = lml.model_real_evaluating(db_real_bitgenia,1,0,0,"no",1,precision=True)
cap_bitgenia_metrics_setpreciso = lml.model_real_evaluating(db_real_bitgenia,0,1,0,"no",1,precision=True)
sim_bitgenia_metrics_setpreciso = lml.model_real_evaluating(db_real_bitgenia,0,0,1,"no",1,precision=True)
j_bitgenia_metrics_setpreciso = lml.model_real_evaluating(db_real_bitgenia,0,0,0,"si",4,precision=True)


# clinvar_metrics_setpreciso = lml.model_real_evaluating(db_real_clinvar,1,1,1,"no",1,precision=True)
##  }}}

##{{{


porcentaje_incluido_bitgenia = [lml.percent_below_x(bitgenia_metrics,i) for i in range(1,2000,50)]
porcentaje_incluido_bitgenia_esp = [lml.percent_below_x(esp_bitgenia_metrics,i) for i in range(1,2000,50)]
porcentaje_incluido_bitgenia_cap = [lml.percent_below_x(cap_bitgenia_metrics,i) for i in range(1,2000,50)]
porcentaje_incluido_bitgenia_sim = [lml.percent_below_x(sim_bitgenia_metrics,i) for i in range(1,2000,50)]
porcentaje_incluido_bitgenia_j = [lml.percent_below_x(j_bitgenia_metrics,i) for i in range(1,2000,50)]

num_hc = range(1,2000,50)

# porcentaje_incluido_clinvar = [lml.percent_below_x(clinvar_metrics,i) for i
        # in range(1,2000,50)]


porcentaje_incluido_bitgenia_setpreciso = [lml.percent_below_x(bitgenia_metrics_setpreciso,i) for i in range(1,2000,50)]
porcentaje_incluido_bitgenia_esp_setpreciso = [lml.percent_below_x(esp_bitgenia_metrics_setpreciso,i) for i in range(1,2000,50)]
porcentaje_incluido_bitgenia_cap_setpreciso = [lml.percent_below_x(cap_bitgenia_metrics_setpreciso,i) for i in range(1,2000,50)]
porcentaje_incluido_bitgenia_sim_setpreciso = [lml.percent_below_x(sim_bitgenia_metrics_setpreciso,i) for i in range(1,2000,50)]
porcentaje_incluido_bitgenia_j_setpreciso = [lml.percent_below_x(j_bitgenia_metrics_setpreciso,i) for i in range(1,2000,50)]



# porcentaje_incluido_clinvar_setpreciso = [lml.percent_below_x(clinvar_metrics_setpreciso,i) for i in range(1,2000,50)]

with plt.style.context(['science','ieee','nature']):
    fig, (ax1,ax2) = plt.subplots(1,2)
    # ax.plot(x,y)
    # ax.set_xlim(0,2000)
    # ax1.plot(num_hc,porcentaje_incluido_clinvar,label='clinvar')
    ax1.plot(num_hc,porcentaje_incluido_bitgenia,label='bitgenia')
    ax1.plot(num_hc,porcentaje_incluido_bitgenia_esp,label='bitgenia $\\frac{j}{j+\\log(k)}$')
    ax1.plot(num_hc,porcentaje_incluido_bitgenia_cap,label='bitgenia $\\frac{j}{j+i}$')
    ax1.plot(num_hc,porcentaje_incluido_bitgenia_sim,label='bitgenia $\\frac{j}{j+i+\\log(k)}$')
    ax1.plot(num_hc,porcentaje_incluido_bitgenia_j,color='cyan',label='bitgenia $j$',alpha=.5)
    ax1.set_xlabel('Top N\%')
    ax1.set_ylabel('Porcentaje capturado')
    ax1.legend(fontsize=4.5)
    ax1.tick_params(axis='both', which='major', labelsize=3)
    ax1.set_title('Set vago')
    ax1.vlines(500,0,0.8,linestyles='dashed',colors='grey')
    ax1.hlines(0.8,0,500,linestyles='dashed',colors='grey')

    # ax2.plot(num_hc,porcentaje_incluido_clinvar_setpreciso,label='clinvar')
    ax2.plot(num_hc,porcentaje_incluido_bitgenia_setpreciso,label='bitgenia')
    ax2.plot(num_hc,porcentaje_incluido_bitgenia_esp_setpreciso,label='bitgenia $\\frac{j}{j+\\log(k)}$')
    ax2.plot(num_hc,porcentaje_incluido_bitgenia_cap_setpreciso,label='bitgenia $\\frac{j}{j+i}$')
    ax2.plot(num_hc,porcentaje_incluido_bitgenia_sim_setpreciso,label='bitgenia $\\frac{j}{j+i+\\log(k)}$')
    ax2.plot(num_hc,porcentaje_incluido_bitgenia_j_setpreciso,color='cyan',label='bitgenia $j$',alpha=0.5)
    ax2.set_xlabel('Top N\%')
    ax2.legend(fontsize=4.5)
    ax2.tick_params(axis='both', which='major', labelsize=3)
    ax2.set_title('Set preciso')
    ax2.vlines(500,0,0.8,linestyles='dashed',colors='grey')
    ax2.hlines(0.8,0,500,linestyles='dashed',colors='grey')


## }}}

## {{{ probando alpha beta y gamma para el heatmap

abg = []

for i in np.arange(0,1.1,0.05):
    for j in np.arange(0,1.1,0.05):
        for k in np.arange(0,1.1,0.05):
            if i+j+k == 1:
                abg.append((i,j,k))

## }}}


## {{{
db = f'{PATH}data/real/bitgenia.json'

with open(db,'r') as f:
    db_real_bitgenia = json.load(f)

# This is a placeholder for your model, replace with your actual model
def evaluate_model(alpha, beta, gamma):
    bitgenia_metrics = lml.model_real_evaluating(db_real_bitgenia,alpha,beta,gamma,"si",1,precision=False)
    percent_below_50 = lml.percent_below_x(bitgenia_metrics,50)
    return percent_below_50

##}}}

## {{{

# Generate range of values for alpha, beta, and gamma
alpha_values = np.linspace(0, 1, 10)
beta_values = np.linspace(0, 1, 10)

# Prepare arrays to hold parameter values and corresponding accuracies
alphas = []
betas = []
gammas = []
accuracies = []

# Calculate model accuracy for each combination of alpha, beta and gamma
for alpha in alpha_values:
    for beta in beta_values:
        # Ensure alpha + beta + gamma <= 1
        if alpha + beta <= 1:
            gamma = 1 - alpha - beta
            accuracy = evaluate_model(alpha, beta, gamma)

            # Store values
            alphas.append(alpha)
            betas.append(beta)
            gammas.append(gamma)
            accuracies.append(accuracy)

# Convert lists to numpy arrays for use with Matplotlib
alphas = np.array(alphas)
betas = np.array(betas)
gammas = np.array(gammas)
accuracies = np.array(accuracies)

## }}}

## {{{
# Create 3D scatter plot
with plt.style.context(['science','ieee','nature']):
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    sc = ax.scatter(alphas, betas, gammas, c=accuracies,
            cmap='RdBu',alpha=1)


    ax.set_xlabel('$\\alpha$', labelpad=-13, size=7)
    ax.set_ylabel('$\\beta$', labelpad=-13, size=7)
    ax.set_zlabel('$\\gamma$', labelpad=-13, size=7)

    ax.tick_params(axis='both', which='major', labelsize=4,
            pad=-5,labelcolor='#4A4A4A')

    # Add a color bar to the plot

    cbar = plt.colorbar(sc)
    cbar.ax.tick_params(labelsize=4)  # Make colorbar labels smaller
    cbar.set_label('Porcentaje en el Top 50',
            size=5,rotation=180+90,labelpad=5.5)  # Label the colorbar

    plt.title('$\\alpha j  - \\beta i -\\gamma k$',size=4.5)

    plt.show()


## }}}





## {{{ corriendo control de métricas y gráfico de accuracy vs top

with open(f'{PATH}data/clinical_cases/synthetic_a_over_a_plus_c.json','r') as f:
    synth_clinical_cases = json.load(f)
    synth_clinical_cases_genes = random.sample([gene for gene in
        synth_clinical_cases],165)
    synth_clinical_cases = {gene:synth_clinical_cases[gene][:5] for gene in
            synth_clinical_cases_genes}

with open(f'{PATH}data/clinical_cases/bitgenia.json','r') as f:
    bitgenia_clinical_cases = json.load(f)

with open(f'{PATH}data/simulated/vague_gene_phenotype_dict.json','r') as f:
    gold_standard_dataset = json.load(f)


# control = lml.model_real_evaluating(synth_clinical_cases,1,0,0,
        # "si",0,precision=False)
ranking_freq_realtiva = lml.model_real_evaluating(bitgenia_clinical_cases,1,0,0,"si",0,precision=False)
# orpha = orpha.v2_model_evaluation(bitgenia_clinical_cases)

## }}}

##{{{
frq_relativa  = [24,
 9,
 932,
 5,
 121,
 38,
 730,
 4,
 394,
 116,
 181,
 23,
 2702,
 5,
 11,
 1,
 751,
 8,
 2238,
 242,
 2041,
 28,
 1034,
 27,
 612,
 75,
 22,
 110,
 71,
 1,
 56,
 350,
 1119,
 1583,
 303,
 6,
 502,
 46,
 3,
 118,
 22,
 6,
 104,
 694,
 1,
 5,
 20,
 40,
 12,
 1,
 6,
 1095,
 148,
 30,
 122,
 84,
 92,
 154,
 318,
 653,
 2,
 135,
 1,
 14,
 356,
 138,
 347,
 390,
 398,
 145,
 1,
 294,
 1,
 67,
 332,
 12,
 1844,
 3,
 17,
 92,
 1,
 39,
 12,
 110,
 512,
 659,
 6,
 10,
 100,
 200,
 15,
 26,
 361,
 96,
 263,
 1717,
 332,
 20,
 2055,
 325,
 344,
 3,
 423,
 384,
 1,
 174,
 326,
 1277,
 70,
 168,
 30,
 20,
 1220,
 147,
 59,
 334,
 273,
 70,
 707,
 6,
 876,
 1295,
 579,
 93,
 59,
 1427,
 14,
 744]

a_over_a_plus_c = [10,41,791,13,5,
 38,
 729,
 5,
 752,
 51,
 86,
 14,
 2770,
 3,
 4,
 4,
 419,
 8,
 1515,
 240,
 311,
 40,
 256,
 58,
 61,
 107,
 41,
 109,
 103,
 1,
 56,
 118,
 586,
 1722,
 303,
 11,
 500,
 508,
 25,
 118,
 7,
 5,
 66,
 693,
 2,
 32,
 27,
 17,
 12,
 1,
 6,
 510,
 61,
 10,
 81,
 103,
 81,
 153,
 335,
 366,
 40,
 135,
 1,
 3,
 356,
 137,
 345,
 390,
 398,
 256,
 1,
 292,
 4,
 3,
 269,
 34,
 96,
 14,
 17,
 5,
 1,
 39,
 42,
 110,
 673,
 490,
 11,
 16,
 160,
 7,
 15,
 2,
 104,
 104,
 1,
 1729,
 181,
 10,
 161,
 314,
 244,
 3,
 63,
 390,
 1,
 395,
 15,
 1350,
 29,
 78,
 31,
 19,
 1356,
 147,
 9,
 301,
 261,
 97,
 169,
 3,
 122,
 1469,
 577,
 243,
 59,
 343,
 4,
 743]

## }}}



## {{{

top_range = [i for i in range(1,101)]
orpha_acc = [lml.percent_below_x(orpha,i) for i in range(1,101)]
control_acc = [lml.percent_below_x(control,i) for i in range(1,101)]
ranking_acc = [lml.percent_below_x(ranking_freq_realtiva,i) for i in range(1,101)]
a_over_a_plus_c_acc = [lml.percent_below_x(a_over_a_plus_c,i) for i in range(1,101)]


## }}}

##{{{
sin_peso = lml.model_real_evaluating(bitgenia_clinical_cases,1,0,0,"si",1,precision=False)
sin_peso_acc = [lml.percent_below_x(sin_peso,i) for i in np.linspace(1,100,100)]
## }}}

## {{{
with plt.style.context(['science','ieee','nature']):
    fig, ax = plt.subplots()
    ax.plot(range,orpha_acc,label='Peso empírico (ORPHA)')
    ax.plot(range,control_acc,label='Control')
    ax.plot(range,ranking_acc,label='Promiscuidad')
    ax.plot(range,a_over_a_plus_c_acc,label='a/(a+c)')
    ax.plot(range,sin_peso_acc,label='Sin peso',linestyle='--')
    ax.set_xlabel('TOP\\#')
    ax.set_ylabel('\\% dentro del TOP\\#')
    ax.legend(fontsize=5)

## }}}




## {{{

"""
Acá voy a hacer los plots para el informe de genphenia, el gráfico sin peso el
gráfico con peso empírico y el gráfico con penalización hpo.

"""

with open(f'{PATH}output/rendimientos_de_orpha_specific_filtered.json','r') as f:
    filtered_rendimientos = json.load(f)





rendimientos = [    list(filtered_rendimientos['sin_pesos'].values()),
                   list(filtered_rendimientos['orpha_sin_penalizacion'].values()),
                   list(filtered_rendimientos['orpha_penalizando'].values())]

recall_at_N = []

for model in rendimientos:
    model_recall = []
    for i in range(1,101):
        recall = lml.percent_below_x(model,i)
        model_recall.append(recall)
    recall_at_N.append(model_recall)

with plt.style.context(['science','ieee','nature']):
    fig,ax = plt.subplots()
    ax.plot(range(1,101),recall_at_N[0],label='Sin pesos')
    ax.plot(range(1,101),recall_at_N[1],label='Peso empírico')
    ax.plot(range(1,101),recall_at_N[2],label='Peso empírico + \nPenalización HPO')
    ax.set_xlabel('TOP\\#')
    ax.set_ylabel('\\% dentro del TOP\\#')
    ax.legend(fontsize=5)




## Lo de abajo crea los boxplots de distribución de rankings
# rendimientos = [[value for value in category if value <= 400] for category in rendimientos]

# Creating the boxplots using the provided template
# with plt.style.context(['science', 'ieee', 'nature']):
    # fig, ax = plt.subplots()
    # ax.boxplot(rendimientos, meanline=True,
               # flierprops=dict(markerfacecolor='r',
                   # markersize=0.4),showfliers=False)
    # ax.set_xlabel('', fontsize=4)
    # ax.set_ylabel('Ranking', fontsize=4)
    # ax.set_title(f'Distribución de rankings para cada feature | Bitgenia', fontsize=4)
    # # plt.subplots_adjust(top=0.9, bottom=0.15,
                        # # left=0.176, right=0.952,
                        # # hspace=0., wspace=0.1)

    # ax.tick_params(axis='both', which='major', labelsize=4)
    # plt.xticks([1, 2, 3], ['Sin pesos', 'Pesos empíricos', 'Pesos empíricos + \n penalización'], rotation=45, fontsize=4)

# # Displaying the plot
# plt.show()




## }}}



## {{{ Distribución de distancias en el grafo HPO de los fenotipos no encontrados a la primera
with open(f'{PATH}data/clinical_cases/bitgenia.json','r') as f:
    clinical_cases = json.load(f)


distances = orpha.v2_model_evaluation(clinical_cases)[1]

## }}}


## {{{ plot de lo de arriba
distances_cleaned = [i for sublist in distances for i in sublist if i>=1]

with plt.style.context(['science','ieee','nature']):
    fig,ax = plt.subplots()
    ax.hist(distances_cleaned,bins=8)
    ax.set_xlabel('Distancia al término más específico')
    ax.set_ylabel('Frecuencia')
    ax.set_title('Fenotipos observados en casos clínicos de Clinvar',fontsize=5)


##}}}


##{{{ Distribución de pesos empíricos en los fenotipos observados
with open(f'{PATH}data/clinical_cases/bitgenia.json','r') as f:
    clinical_cases = json.load(f)

weights = orpha.v2_model_evaluation(clinical_cases)[2]

## }}}

##{{{ Distribución de pesos empíricos en los fenotipos observados
with open(f'{PATH}data/clinical_cases/clinvar.json','r') as f:
    clinical_cases_clinvar = json.load(f)

weights_clinvar = orpha.v2_model_evaluation(clinical_cases_clinvar)[2]

## }}}

##{{{
weight_distances_cleaned = [i for sublist in weights_clinvar for i in sublist if i>=1]

with plt.style.context(['science','ieee','nature']):
    fig,ax = plt.subplots()
    ax.hist(weight_distances_cleaned)
    ax.set_xlabel('Pesos empíricos de ORPHA')
    ax.set_ylabel('Frecuencia')
    ax.set_title('Pesos empíricos en fenotipos observados de Clinvar',fontsize=5)
    ax.set_xticks([2,17,54,89])
    tick_labels = [2.5, 17.0, 54.5, 89.5]
    ax.set_xticklabels(tick_labels)




## }}}


##{{{
distances_cleaned = [i for sublist in weights for i in sublist if i>=1]

with plt.style.context(['science','ieee','nature']):
    fig,ax = plt.subplots()
    ax.hist(distances_cleaned,bins=6)
    ax.set_xlabel('Pesos empíricos de ORPHA')
    ax.set_ylabel('Frecuencia')
    ax.set_title('Pesos empíricos en fenotipos observados de Bitgenia',fontsize=5)
    ax.set_xticks([2,17,54,89])
    tick_labels = [2.5, 17.0, 54.5, 89.5]
    ax.set_xticklabels(tick_labels)





## }}}
