# Here are the main codes for running scGMM-VGAE
# These main codes are used for clustering the labeled data
# The codes here are further modified based on Eric Lin's codes (https://github.com/ericlin1230/scGMM-VGAE)
# Lin., Liu., et al., 2022 (In preparation)

import numpy as np
import torch
import scipy.sparse as sp
from model_CUDA_labeled import GMM_VGAE
from preprocessing import load_data, sparse_to_tuple, preprocess_graph
import time

# Dataset Name
dataset = "Zhalf1"
print("Zhalf1")
nClusters = 4
adj, features, labels = load_data('Zhalf1', './data/Zhalf1', True)
features_new = features.toarray()
num_nodes = features.shape[0]
num_features = features.shape[1]

# Network parameters
num_neurons = 24
embedding_size = 12
save_path = "./results/"

# Pretraining parameters
epochs_pretrain = 200
lr_pretrain = 0.01

# Clustering parameters
epochs_cluster = 200
lr_cluster = 0.01

# Configure the device to cuda
torch.set_default_tensor_type('torch.cuda.FloatTensor') 
device = torch.device("cuda")

# Data processing
adj = adj - sp.dia_matrix((adj.diagonal()[np.newaxis, :], [0]), shape=adj.shape)
adj.eliminate_zeros()
adj_norm = preprocess_graph(adj)
features = sparse_to_tuple(features.tocoo())
num_features = features[2][1]
pos_weight_orig = float(adj.shape[0] * adj.shape[0] - adj.sum()) / adj.sum()
norm = adj.shape[0] * adj.shape[0] / float((adj.shape[0] * adj.shape[0] - adj.sum()) * 2)
adj_label = adj + sp.eye(adj.shape[0])
adj_label = sparse_to_tuple(adj_label)
adj_norm = torch.cuda.sparse.FloatTensor(torch.LongTensor(adj_norm[0].T).cuda(), torch.FloatTensor(adj_norm[1]).cuda(), torch.Size(adj_norm[2]))
adj_label = torch.cuda.sparse.FloatTensor(torch.LongTensor(adj_label[0].T).cuda(), torch.FloatTensor(adj_label[1]).cuda(), torch.Size(adj_label[2]))
features = torch.cuda.sparse.FloatTensor(torch.LongTensor(features[0].T).cuda(), torch.FloatTensor(features[1]).cuda(), torch.Size(features[2]))
weight_mask_orig = adj_label.to_dense().view(-1) == 1
weight_tensor_orig = torch.ones(weight_mask_orig.size(0))
weight_tensor_orig[weight_mask_orig] = pos_weight_orig

print("start")
# Start the timer
start = time.perf_counter()
# Training
acc_array = []

from random import randint
ress = []
seed = 1 # Change seed here
network = GMM_VGAE(adj = adj_norm , num_neurons=num_neurons, num_features=num_features, embedding_size=embedding_size, nClusters=nClusters, activation="Sigmoid", seed=seed)
network.to(device)
acc_array = network.pretrain(adj_norm, features, adj_label, labels, weight_tensor_orig, norm, optimizer="Adam", epochs=epochs_pretrain, lr=lr_pretrain, save_path=save_path, dataset=dataset)
network.to(device)
res, y_pred, y = network.train([], adj_norm, features, adj_label, labels, weight_tensor_orig, norm, optimizer="Adam", epochs=epochs_cluster, lr=lr_cluster, save_path=save_path, dataset=dataset, features_new=features_new)
end = time.perf_counter()

print(f"Total time: {end - start:0.4f} seconds")

print(seed)

