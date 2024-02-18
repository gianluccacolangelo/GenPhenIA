"""
Acá voy a ver si orpha contiene a omim
"""

## {{{
import pandas as pd
PATH = "/home/brainy/Desktop/1ercuatri2023/Tesis/GenPhenIA/"
import xml.etree.ElementTree as ET
import csv
## }}}



##{{{
tree = ET.parse(f'{PATH}data/ORPHA/alignments_omim_orpha.xml')
root = tree.getroot()

def print_structure(element, level=0):
    print('  ' * level + element.tag)
    i = 0
    for child in element:
        i+=1
        print_structure(child, level + 1)
        if i==3:
            break

## }}


## {{{
mim_orpha_dict = {} # key: mim, value: orpha
counter=0
for disorder in root.findall('.//Disorder'):
    ext_ref_list = disorder.find('ExternalReferenceList')
    orphacode = disorder.find('OrphaCode').text

    if ext_ref_list is not None:
        # Iterate over all 'ExternalReference' elements in the list
        for ext_ref in ext_ref_list.findall('ExternalReference'):
            # Extract and print information from each 'ExternalReference'
            source = ext_ref.find('Source').text if ext_ref.find('Source') is not None else None
            reference = ext_ref.find('Reference').text if ext_ref.find('Reference') is not None else None
            # print(f'Source: {source}, Reference: {reference}')
            # print(orphacode)
            if source=='OMIM':
                counter+=1
                mim_orpha_dict[reference] = [orphacode]
## }}}


## {{{




# Path to your file
file_path = f'{PATH}data/phenotype.hpoa'

# List to store MIM codes
mim_codes = []

# Open the file and create a csv reader


with open(file_path, 'r') as file:
    reader = csv.DictReader(file, delimiter='\t')

    # Iterate over each row in the file
    for row in reader:
        # The MIM code is in the '#description: "HPO annotations for rare diseases [8120: OMIM; 47: DECIPHER; 4264 ORPHANET]"' column
        database_id = row.get('#description: "HPO annotations for rare diseases [8120: OMIM; 47: DECIPHER; 4264 ORPHANET]"', '')

        # Check if the ID is a MIM code
        if database_id.startswith('OMIM:'):
            # Extract the MIM code and add it to the list
            mim_code = database_id.replace('OMIM:', '')
            mim_codes.append(mim_code)


## }}}


## {{{
total_mim_codes = list(set(mim_codes))

total_in = 0
total_out = 0
out_mim = []
in_mim = []
for mim in total_mim_codes:
    if mim in mim_orpha_dict.keys():
        total_in+=1
        in_mim.append(mim)
    else:
        total_out+=1
        out_mim.append(mim)
print(total_in/len(total_mim_codes))


## }}}


## {{{


phen_disease_tree = ET.parse(f'{PATH}data/ORPHA/phen_disease.xml')
phen_disease_root = phen_disease_tree.getroot()

print_structure(phen_disease_root)

for hpo_disorder in phen_disease_root.findall('.//HPODisorderSetStatus'):
    disorder = hpo_disorder.find('Disorder')


## }}}
