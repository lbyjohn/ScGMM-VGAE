# Main codes of Seurat clustering are in Seurat_Clustering.R
# Following codes are to evaluate the clustering performance of Seurat

from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import normalized_mutual_info_score
from sklearn.metrics import silhouette_score

# Load the Seurat clustering results
import pyreadr
topfeat = pyreadr.read_r('DSim/DSim_10.rds')
res = topfeat[None].iloc[:, 0]
res = list(res)

# Load the truth labels
import pickle as pkl
file = open('DSim/DSim_label.pkl', 'rb')
realres = pkl.load(file)

# Calculate ARI, NMI, and SC values
ARI=adjusted_rand_score(realres, res)
print(ARI)
NMI=normalized_mutual_info_score(realres, res)
print(NMI)
import pandas as pd
df2 = pd.read_csv (r'DSim/DSim_10.csv')
df2 = df2.iloc[: ,1:]
df2 = df2.transpose()
features_new = df2.values
sc = silhouette_score(features_new, res, metric="cosine")
print(sc)