"""
En este script modificamos los sets reales para recuperar a los términos hpo en
su forma original "HP:XXXXXXXX" y no solo como números. Para que sea
consistente con los sets simulados y pueda ser utilizado con las funciones
ya definidas
"""

## {{{ IMPORTACIONES
import json
## }}}



## {{{


import json

def format_hpo_id(hpo_id):
    # format the HPO id to have 'HP:' prefix and seven digits
    return f'HP:{str(hpo_id).zfill(7)}'

def main():
    with open('clinvar.json') as f:
        data = json.load(f)

    for gene_id, hpo_ids in data.items():
        data[gene_id] = [format_hpo_id(hpo_id) for hpo_id in hpo_ids]

    with open('output.json', 'w') as f:
        json.dump(data, f)

if __name__ == "__main__":
    main()
##  }}}
