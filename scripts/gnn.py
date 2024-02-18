###{{{ importaciones
import dgl
import torch
import torch.nn as nn
import torch.nn.functional as F

###}}}


##{{{
# Load the Zachary's Karate Club dataset
def load_karate_club():
    dataset = dgl.data.KarateClubDataset()
    graph = dataset[0]
    # Add self-loop to each node (required for GCN)
    graph = dgl.add_self_loop(graph)
    # Prepare labels and masks for semi-supervised learning
    labels = graph.ndata['label']
    mask = torch.zeros(len(labels), dtype=torch.bool)
    mask[[0, 33]] = True  # Only the instructor and the club president nodes are labeled
    return graph, labels, mask

graph, labels, train_mask = load_karate_club()
##}}}
