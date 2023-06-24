## {{{ importaciones
from owlready2 import *
PATH = "/home/brainy/Desktop/1ercuatri2023/Tesis/GenPhenIA/"
## }}}


## {{{
my_world = World()
ontology = my_world.get_ontology(f"{PATH}data/ORPHA/hoom_orphanet.owl").load()

## {{{


# Access the 'hoom_orphanet.Association' class

association_class = ontology.search_one(iri="*Association")

associations = list(association_class.subclasses())


for association in associations[:15]:  # Adjust this to explore more or less subclasses
    print(f"Association: {association}")

    # Get the properties and their values
    print(f"association_has_subject: {association.association_has_subject}")
    print(f"association_has_object: {association.association_has_object}")
    print(f"has_frequency: {association.has_frequency}")
    print(f"has_DC_attribute: {association.has_DC_attribute}")
    print(f"has_provenance: {association.has_provenance}")

    # Explore 'Provenance' and its property 'has_evidence'
    if association.has_provenance:
        for provenance in association.has_provenance:
            print(f"Provenance: {provenance}")
            try:
                print(f"has_evidence: {provenance.has_evidence}")
            except:
                continue

    print()  # Empty line for better readability

## }}}
