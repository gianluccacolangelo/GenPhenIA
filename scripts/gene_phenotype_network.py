"""
En este script vamos a construir y analizar una red gen-fenotipo utilizando datos de Hu-
man Phenotype Ontology (HPO) y las bases de datos de OMIM, con el objetivo
de identificar comunidades de genes asociados con fenotipos similares
"""


## {{{
import networkx as nx
import json
import community as community_louvain
import numpy as np
PATH = "/home/brainy/Desktop/1ercuatri2023/Tesis/GenPhenIA/"
# Load your gene-phenotype associations
with open(f'{PATH}data/simulated/gene_phenotype_dict.json', 'r') as f:
    gene_to_phenotype = json.load(f)
## }}}

## {{{

# Create a new bipartite graph
B = nx.Graph()

# Add nodes with the bipartite attribute
for gene in gene_to_phenotype.keys():
    B.add_node(gene, bipartite=0)
for phenotype in set([p for phenotypes in gene_to_phenotype.values() for p in phenotypes]):
    B.add_node(phenotype, bipartite=1)

# Add edges
for gene, phenotypes in gene_to_phenotype.items():
    for phenotype in phenotypes:
        B.add_edge(gene, phenotype)

# Check if the graph is connected
print(nx.is_connected(B))

# Get the set of genes and phenotypes
genes = set(n for n, d in B.nodes(data=True) if d['bipartite']==0)
phenotypes = set(B) - genes



# Create the weighted projected graph
G = nx.bipartite.weighted_projected_graph(B, genes)

# Apply a threshold
# threshold = 30  # adjust this value based on your specific needs
# edges_below_threshold = [(u, v) for u, v, d in G.edges(data=True) if d['weight'] < threshold]
# G.remove_edges_from(edges_below_threshold)

# Print some statistics about the graph
# print(nx.info(G))

## }}}

##{{{
largest_component = max(nx.connected_components(G), key=len)
components_sorted_by_size = sorted(nx.connected_components(G), key=len, reverse=True)
## }}}

## {{{

def compute_metrics(G, threshold):
    # Remove edges with weight less than the threshold
    G_thresholded = G.copy()
    for u, v, data in G.edges(data=True):
        if data['weight'] < threshold:
            G_thresholded.remove_edge(u, v)

    # Find the largest component
    largest_component = max(nx.connected_components(G_thresholded), key=len)

    # Compute the modularity
    partition = community_louvain.best_partition(G_thresholded)
    modularity = community_louvain.modularity(partition, G_thresholded)

    return modularity, len(largest_component)

def find_optimal_threshold(G, thresholds):
    # Compute the metrics for each threshold
    metrics = [compute_metrics(G, t) for t in thresholds]

    # Find the threshold that maximizes the modularity and the size of the largest component
    optimal_threshold = thresholds[np.argmax([m[0] + m[1] for m in metrics])]

    return optimal_threshold

## }}}


## {{{
thresholds = np.linspace(2, 10, 9)  # adjust this range based on your specific needs
optimal_threshold = find_optimal_threshold(G, thresholds)
print('Optimal threshold:', optimal_threshold)
##}}}




## {{{

largest_cc = max(nx.connected_components(G), key=len)

# Create a subgraph of G consisting only of this component
Gc = G.subgraph(largest_cc)

# Compute the best partition of this subgraph
partition = community_louvain.best_partition(Gc)
modularity = community_louvain.modularity(partition, Gc)
print(modularity)
# Create a color map (one color for each community)
colors = [partition[node] for node in Gc.nodes]

# Draw the nodes with a smaller size
nx.draw_networkx_nodes(Gc, pos=nx.spring_layout(Gc), node_color=colors,
        node_size=40,
        alpha=.5)

# Draw the edges
nx.draw_networkx_edges(Gc, pos=nx.spring_layout(Gc), alpha=0.01,width=1)

plt.show()

## }}}

## {{{
#Definamos una función que, dada una red y una partición, nos devuelva un plot del clustering de cada nodo
def silhouette(Red,particion): # La función nos pide la Red y la partición como diccionario
  S=[]
  limites=[0]
  # Recorro los clusters sin repetir
  for cluster in set(particion.values()):
    #Filtro los nodos que pertenecen a este cluster
    nodos_en_cluster = [nodo for (nodo, value) in particion.items() if value == cluster]
    S_cluster=[]
    # Recorro los nodos del cluster
    for nodo in nodos_en_cluster:
      distancias_dentro=[]
      distancias_fuera=[]
      # Recorro los nodos del mismo cluster
      for nodo_en_cluster in nodos_en_cluster:
        if nodo != nodo_en_cluster:
          # Calculo y guardo la distancia, si no es consigo mismo
          distancias_dentro.append(nx.shortest_path_length(Red, source=nodo, target=nodo_en_cluster, weight=None))
      # Recorro los nodos de los otros clusters
      for nodo_fuera in Red.nodes():
        if particion[nodo_fuera] != cluster:
          # Calculo y guardo la distancia
          distancias_fuera.append(nx.shortest_path_length(Red, source=nodo, target=nodo_fuera, weight=None))
      # Calculo la distancia media para los del mismo cluster
      distancia_media_dentro = np.mean(distancias_dentro)
      # Calculo la distancia mínima para los nodos fuera del cluster
      distancia_mean_fuera = np.mean(distancias_fuera)
      # Calculo y guardo la Silhouette del nodo
      S_cluster.append((distancia_mean_fuera - distancia_media_dentro)/np.max([distancia_mean_fuera,distancia_media_dentro]))
    # Ordeno las Silhouette del mismo cluster por valor, para graficar lindo
    S_cluster=sorted(S_cluster)
    # Me guardo en qué nodo termina cada cluster, para graficar clusters por colores
    limites.append(len(S_cluster)+limites[-1])
    # Agrego las Silhouette de este cluster a la lista de todas
    S = S + S_cluster
  # Calculo la Silhouette media
  S_media = np.mean(S)
  # Grafico todas con colores por clusters
  for i in range(len(limites)-1):
    plt.barh(range(limites[i],limites[i+1]), S[limites[i]:limites[i+1]])
#    plt.plot(range(limites[i],limites[i+1]),S[limites[i]:limites[i+1]])
  plt.axvline(S_media)
  return S, S_media

## }}}

## {{{
s, s_media = silhouette(Gc, partition)
plt.title(nombres[3])
plt.show()
## }}}
