# The main code CellVGAE_Clustering is a shell script, which is based on the Cell-VGAE algorithm (Buterez et al., 2021)
# The source code of Cell-VGAE can be found at https://github.com/davidbuterez/CellVGAE
# Following codes are to evaluate the clustering performance of Cell-VGAE

# Load the Cell-VGAE clustering results
import numpy as np
clusters = np.load('./DSim_10/hdbscan_clusters/hdbscan-clusters-min_cluster_size=100-min_samples=50.npy')

# Load the truth labels
import pickle as pkl
file = open('DSim_label.pkl', 'rb')
realres = pkl.load(file)

# Calculate ARI
import sklearn.metrics
ARI = sklearn.metrics.adjusted_rand_score(realres, clusters)
print(ARI)

# Calculate NMI
NMI = sklearn.metrics.normalized_mutual_info_score(realres, clusters)
print(NMI)

# Calculate SC
import pandas as pd
df2 = pd.read_csv (r'Rawdata/DSim_10.csv')
df2 = df2.iloc[: ,1:]
df2 = df2.transpose()
features_new = df2.values
sc = sklearn.metrics.silhouette_score(features_new, clusters, metric="cosine")
print(sc)