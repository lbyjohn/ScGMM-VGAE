# This code is based on the DCA + Leiden algorithm (Eraslan et al., 2019)
# Source code of DCA + Leiden can be found at https://github.com/theislab/dca

# Load the data
import pandas as pd
df = pd.read_csv (r'DSim/DSim_10.csv')
df = df.iloc[:,1:]
df = df.transpose()
features_new = df.values

# Pre-processing
import scanpy as sc
adata = sc.AnnData(df.values)
adata

# DCA denoising
from dca.api import dca
dca(adata, threads=10)

# Leiden clustering
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)
sc.tl.tsne(adata)

sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.leiden(adata)

# Load the truth labels
import pickle as pkl
file = open('DSim/DSim_label.pkl', 'rb')
realres = pkl.load(file)

# Performance evaluation (ARI, NMI, SC metrics)
res = list(adata.obs.leiden)
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import normalized_mutual_info_score
ARI=adjusted_rand_score(realres, res)
NMI=normalized_mutual_info_score(realres, res)
print(ARI)
print(NMI)
from sklearn.metrics import silhouette_score
sc = silhouette_score(features_new, res, metric="cosine")
print(sc)
