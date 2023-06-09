"""
En este archivo están los plots de todos los gráficos.
"""

## {{{ importaciones
import sys
sys.path.insert(0,'/home/brainy/Desktop/1ercuatri2023/Tesis/GenPhenIA/src')
import phen_gen_weight_functions as pgw
import linear_model_lab as lml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import json
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

clinvar_metrics = lml.model_real_evaluating(db_real_clinvar,1,1,1,"no",1,precision=False)

bitgenia_metrics_setpreciso = lml.model_real_evaluating(db_real_bitgenia,1,1,1,"no",1,precision=True)

clinvar_metrics_setpreciso = lml.model_real_evaluating(db_real_clinvar,1,1,1,"no",1,precision=True)
##  }}}

##{{{


porcentaje_incluido_bitgenia = [lml.percent_below_x(bitgenia_metrics,i) for i in range(1,2000,50)]
num_hc = range(1,2000,50)

porcentaje_incluido_clinvar = [lml.percent_below_x(clinvar_metrics,i) for i
        in range(1,2000,50)]


porcentaje_incluido_bitgenia_setpreciso = [lml.percent_below_x(bitgenia_metrics_setpreciso,i) for i in range(1,2000,50)]

porcentaje_incluido_clinvar_setpreciso = [lml.percent_below_x(clinvar_metrics_setpreciso,i) for i in range(1,2000,50)]

with plt.style.context(['science','ieee','nature']):
    fig, (ax1,ax2) = plt.subplots(1,2)
    # ax.plot(x,y)
    # ax.set_xlim(0,2000)
    ax1.plot(num_hc,porcentaje_incluido_clinvar,label='clinvar')
    ax1.plot(num_hc,porcentaje_incluido_bitgenia,label='bitgenia')
    ax1.set_xlabel('Top N\%')
    ax1.set_ylabel('Porcentaje capturado')
    ax1.legend()
    ax1.tick_params(axis='both', which='major', labelsize=3)
    ax1.set_title('Set vago')
    ax1.vlines(500,0,0.8,linestyles='dashed',colors='grey')
    ax1.hlines(0.8,0,500,linestyles='dashed',colors='grey')

    ax2.plot(num_hc,porcentaje_incluido_clinvar_setpreciso,label='clinvar')
    ax2.plot(num_hc,porcentaje_incluido_bitgenia_setpreciso,label='bitgenia')
    ax2.set_xlabel('Top N\%')
    ax2.legend()
    ax2.tick_params(axis='both', which='major', labelsize=3)
    ax2.set_title('Set preciso')
    ax2.vlines(500,0,0.8,linestyles='dashed',colors='grey')
    ax2.hlines(0.8,0,500,linestyles='dashed',colors='grey')


## }}}

