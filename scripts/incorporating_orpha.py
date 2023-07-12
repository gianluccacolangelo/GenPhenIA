"""
En este script voy a ir incorporando orpha al modelo de predicción, las
enfermedades, las clasificaciones y los pesos empiricos que brinda.

"""

## {{{ IMPORTACIONES
import pandas as pd
import numpy as np
import xml.etree.ElementTree as ET
import pandas as pd
# sys.path.insert(0,'/home/brainy/Desktop/1ercuatri2023/Tesis/GenPhenIA/')
PATH = "/home/brainy/Desktop/1ercuatri2023/Tesis/GenPhenIA/"
import json
## }}}



## {{{
# Parse the XML file
tree_phen = ET.parse(f'{PATH}/data/ORPHA/phen_disease.xml')
phen_root= tree_phen.getroot()

tree_gene = ET.parse(f'{PATH}/data/ORPHA/gene_disease.xml')
gene_root= tree_gene.getroot()

## }}}



## {{{
with open(f'{PATH}data/OMIM/genemap2.txt', 'r') as f:
    for line in f:
        if line.startswith('# Chromosome'):
            columns = line.strip('#\n').split('\t')
            break

# Load the file
genemap2_df = pd.read_csv(f'{PATH}data/OMIM/genemap2.txt', comment='#', delimiter='\t', names=columns, skiprows=4)

# Drop rows with missing values in 'Gene Symbols' or 'Entrez Gene ID'
genemap2_df = genemap2_df.dropna(subset=['Approved Gene Symbol', 'Entrez Gene ID'])

# Convert 'Entrez Gene ID' to integer
genemap2_df['Entrez Gene ID'] = genemap2_df['Entrez Gene ID'].astype(int)

# Create the dictionaries
gene_to_entrez = pd.Series(genemap2_df['Entrez Gene ID'].values,index=genemap2_df['Approved Gene Symbol']).to_dict()
entrez_to_gene = pd.Series(genemap2_df['Gene Symbols'].values,index=genemap2_df['Entrez Gene ID']).to_dict()

def translate_gene_to_entrez(gene, gene_to_entrez):
    """
    Translate a gene symbol to an Entrez Gene ID.
    """
    return gene_to_entrez.get(gene, None)

def translate_entrez_to_gene(entrez, entrez_to_gene):
    """
    Translate an Entrez Gene ID to a gene symbol.
    """
    return entrez_to_gene.get(entrez, None)
## }}}


## {{{

def frequency_to_weight(frequency):
    """
    Convert a frequency name to a weight.
    """
    if "Always present" in frequency:
        return 100
    elif "Very frequent" in frequency:
        return (99 + 80) / 2
    elif "Frequent" in frequency:
        return (79 + 30) / 2
    elif "Occasional" in frequency:
        return (29 + 5) / 2
    elif "Very rare" in frequency:
        return (4 + 1) / 2
    elif "Excluded" in frequency:
        return 0
    else:
        return None

def parse_xml_to_dict(root,disorder_name=False):
    """
    Parse the XML tree and return a dictionary where each key is a disease,
    and each value is another dictionary where each key is a phenotype and the value is the empirical weight of that phenotype.
    """
    data = {}

    # Iterate over all disorders
    for disorder_set_status in root.find('HPODisorderSetStatusList'):
        disorder = disorder_set_status.find('Disorder')
        if disorder is not None:
            # Get the disease name
            if disorder_name == True:
                disease_name = disorder.find('Name')
                if disease_name is not None:
                    disease_name = disease_name.text

            elif disorder_name == False:
                disease_code = disorder.find('OrphaCode')
                if disease_code is not None:
                    disease_name = disease_code.text

                # Initialize the inner dictionary for this disease
                data[disease_name] = {}

                # Iterate over all HPODisorderAssociation elements (phenotype-frequency pairs) for this disease
                for hpo_disorder_association in disorder.find('HPODisorderAssociationList'):
                    hpo = hpo_disorder_association.find('HPO')
                    hpo_frequency = hpo_disorder_association.find('HPOFrequency')

                    if hpo is not None and hpo_frequency is not None:
                        # Get the phenotype HPO term
                        phenotype = hpo.find('HPOId')
                        if phenotype is not None:
                            phenotype = phenotype.text

                        # Get the frequency and convert it to a weight
                        frequency = hpo_frequency.find('Name')
                        if frequency is not None:
                            weight = frequency_to_weight(frequency.text)

                        # Add the phenotype and its weight to the inner dictionary for this disease
                        if phenotype is not None and weight is not None:
                            data[disease_name][phenotype] = weight

    return data



def parse_xml_to_gene_dict(root,disorder_name=False):

    """
    Parse the XML tree and return a dictionary where each key is a disease,
    and each value is a list of genes associated with that disease.
    """

    data = {}

    # Iterate over all disorders
    for disorder in root.find('DisorderList'):
        # Get the disease name
        if disorder_name == True:
            disease_name = disorder.find('Name')
            if disease_name is not None:
                disease_name = disease_name.text
        elif disorder_name == False:
            disease_name = disorder.find('OrphaCode')
            if disease_name is not None:
                disease_name = disease_name.text

            # Initialize the list for this disease
            data[disease_name] = []

            # Iterate over all DisorderGeneAssociation elements for this disease
            for disorder_gene_association in disorder.find('DisorderGeneAssociationList'):
                gene = disorder_gene_association.find('Gene')

                if gene is not None:
                    # Get the gene symbol
                    gene_symbol = gene.find('Symbol')
                    if gene_symbol is not None:
                        gene_symbol = gene_symbol.text

                    # Add the gene symbol to the list for this disease
                    if gene_symbol is not None:
                        data[disease_name].append(gene_symbol)

    return data


def rank_diseases(observed_phenotypes, phen_disease_dict, threshold=3):
    """
    Rank diseases based on observed phenotypes. From the dictionary where
    diseases are keys and the values are dictionaries of phenotypes and
    weights. The threshold is the minimum number of observed phenotypes that
    must be present in a disease to be considered.
    """
    # Initialize a dictionary to store the total weight for each disease
    total_weights = {}

    # Iterate over all diseases in the phen_disease_dict
    for disease, phenotypes in phen_disease_dict.items():
        # Initialize the total weight for this disease
        total_weight = 0

        # Initialize the count of observed phenotypes for this disease
        count = 0

        # Iterate over all observed phenotypes
        for observed_phenotype in observed_phenotypes:
            # If this phenotype is present in the disease, add its weight to the total weight for this disease
            if observed_phenotype in phenotypes:
                total_weight += phenotypes[observed_phenotype]
                count += 1

        # If the total weight for this disease is nonzero and the count of observed phenotypes is above the threshold,
        # add it to the total_weights dictionary
        if total_weight > 0 and count >= threshold:
            total_weights[disease] = total_weight

    # Sort the diseases by total weight in descending order and return the sorted list
    ranked_diseases = sorted(total_weights.items(), key=lambda item: item[1], reverse=True)
    return ranked_diseases



## }}}



## {{{

# Parse the XML file to a dictionary
phen_disease_dict = parse_xml_to_dict(phen_root,disorder_name=False)

# Print the first few entries in the dictionary to verify that it's correct
for disease, phenotypes in list(phen_disease_dict.items())[:1]:
    print(disease, phenotypes)




observed_phenotypes = ["HP:0011451",
        "HP:0000252",
        "HP:0001298",
        "HP:0002650",
        "HP:0001290",
        "HP:0007360",
        "HP:0001263",
        "HP:0001344",
        "HP:0000316",
        "HP:0000490",
        "HP:0000490",
        "HP:0003196",
        "HP:0009933",
        "HP:0000303",
        "HP:0040080",
        "HP:0002213",
        "HP:0100874",
        "HP:0001167",
        "HP:0040196",
        "HP:0001250",
        "HP:0200134"]

# Parse the XML file to a dictionary
gene_disease_dict = parse_xml_to_gene_dict(gene_root)

# Print the first few entries in the dictionary to verify that it's correct
for disease, genes in list(gene_disease_dict.items())[10:15]:
    print(disease, genes)



ranked_diseases = rank_diseases(observed_phenotypes, phen_disease_dict, threshold=1)
for disease, weight in ranked_diseases[:5]:
    print(f"Disease(OrphaCode): {disease}, Total Weight: {weight}")
    # print(gene_disease_dict[disease])



## }}}

## {{
genes_symbols_list = []
for disease in ranked_diseases:
    try:
        genes_symbols_list.append(gene_disease_dict[disease[0]])
    except:
        continue


## }}





## {{{
with open(f'{PATH}data/clinical_cases/bitgenia.json','r') as f:
    clinical_cases = json.load(f)



def v2_model_evaluation(clinical_cases,top_diseases=50,top_genes=50):
    encontrados = 0
    no_encontrados = 0

    for observed_genes in clinical_cases.keys():
        observed_phenotypes = clinical_cases[observed_genes]
        ranked_diseases = rank_diseases(observed_phenotypes, phen_disease_dict, threshold=1)
        candidate_genes = []

        for disease in ranked_diseases:

            try:
                genes = gene_disease_dict[disease[0]]
            except:
                continue
            if len(genes)>1:
                for gene in genes:
                    gene = translate_gene_to_entrez(gene,gene_to_entrez)
                    if gene not in candidate_genes:
                        candidate_genes.append(str(gene))

            elif len(genes)==1:

                gene = genes[0]
                gene = translate_gene_to_entrez(gene,gene_to_entrez)

                if gene not in candidate_genes:
                    candidate_genes.append(str(gene))

        if observed_genes in candidate_genes:
            encontrados += 1
        else:
            no_encontrados += 1
        try:
            print(f'ranking: {candidate_genes.index(observed_genes)}')
        except:
            print('mm')
    return encontrados,no_encontrados


## }}}
