# This code is based on the scVI algorithm (Lopez et al., 2018)
# Source code of DCA + Leiden can be found at https://github.com/scverse/scvi-tools

import scvi
import scanpy as sc

# Load the data
df = scvi.data.read_csv('DSim/DSim_10.csv')
df = df.transpose()

# Train the scVI model
scvi.model.SCVI.setup_anndata(df)
model = scvi.model.SCVI(df)
model
model.train()
latent = model.get_latent_representation()
df.obsm["X_scVI"] = latent
df.layers["scvi_normalized"] = model.get_normalized_expression(
    library_size=10e4
)

# Clustering
sc.tl.pca(df, svd_solver='arpack')
sc.pp.neighbors(df, n_pcs=30, n_neighbors=20)
sc.tl.umap(df, min_dist=0.3)
df.obs
sc.tl.leiden(df, key_added="leiden_scVI", resolution=0.5)
sc.pl.umap(
    df,
    color=["leiden_scVI"],
    frameon=False,
)
df.obs.leiden_scVI

# Load the truth label
import pickle as pkl
file = open('DSim/DSim_label.pkl', 'rb')
realres = pkl.load(file)

# Performance evaluation (ARI, NMI, SC metrics)
import sklearn.metrics
ARI=sklearn.metrics.adjusted_rand_score(realres, list(df.obs.leiden_scVI))
NMI=sklearn.metrics.normalized_mutual_info_score(realres, list(df.obs.leiden_scVI))
print(ARI)
print(NMI)
import pandas as pd
df2 = pd.read_csv (r'DSim/DSim_10.csv')
df2 = df2.iloc[: ,1:]
df2 = df2.transpose()
features_new = df2.values
sc = sklearn.metrics.silhouette_score(features_new, list(df.obs.leiden_scVI), metric="cosine")
print(sc)


