## {{{
import json
from collections import defaultdict
## }}}

## {{{
def main():
    # Load data from JSON file
    with open('phenotype_gene_dict.json') as f:
        data = json.load(f)

    # Invert dictionary
    inverted_data = defaultdict(list)
    for hpo_term, genes in data.items():
        for gene in genes:
            inverted_data[gene].append(hpo_term)

    # Write the inverted data to a new JSON file
    with open('vague_gene_phenotype_dict.json', 'w') as f:
        json.dump(inverted_data, f)

if __name__ == "__main__":
    main()
## }}}
