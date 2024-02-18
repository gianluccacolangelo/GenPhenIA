## {{{ importaciones
from owlready2 import *
PATH = "/home/brainy/Desktop/1ercuatri2023/Tesis/GenPhenIA/"
## }}}


## {{{
my_world = World()
hoom = my_world.get_ontology(f"{PATH}data/ORPHA/hoom_orphanet.owl").load()
## }}}

## {{{


# Access the 'hoom_orphanet.Association' class

hoom_association_class = hoom.search_one(iri="*Association")


associations = list(hoom_association_class.subclasses())


for association in associations[:15]:  # Adjust this to explore more or less subclasses
    print(f"Association: {association}")

    # Get the properties and their values
    print(f"association_has_subject: {association.association_has_subject}")
    print(f"association_has_object: {association.association_has_object}")
    subj = association.association_has_subject
    print(dir(subj[0]))
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


## {{{

for cls in hoom.classes():
    if 'Orpha:1934' in str(cls):
        print(cls)

for prop in hoom.properties():
    if '1934' in str(prop):
        print(prop)


for cls in ontology.classes():
    for label in cls.label:
        if 'OMIM' in str(label):
            print(f"Class: {cls}, Label: {label}")
    for comment in cls.comment:
        if 'OMIM' in str(comment):
            print(f"Class: {cls}, Comment: {comment}")

## }}}
