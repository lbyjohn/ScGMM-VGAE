# The following codes are to generate the .allx, .ally, and .graph files
# Those files are the inputs for scGMM-VGAE
# .allx file is the transformed processed feature matrix
# .ally file is the transformed truth label
# .graph file is the transformed cell-cell adjacency matrix

# Load the processed feature matrix and truth labels
import pandas as pd
from scipy.sparse import coo_matrix
df = pd.read_csv (r'B3Sim/B3Sim_1.csv')
df2 = pd.read_csv (r'B3Sim/B3Sim_label.csv')
import numpy as np
import scipy
import pickle as pkl

# Transform feature matrix into .allx file
dfv = df
dfv = dfv.iloc[: ,1:]
dff = scipy.sparse.csr_matrix(dfv.values)
dff = dff.transpose()
pkl.dump(dff, open("B3Sim/ind.B3Sim1.allx", "wb"))

# Transform truth labels into .ally file
unique1 = list(set(df2['assigned_cluster']))
dict1 = {}
for i in range(len(unique1)):
    dict1[unique1[i]] = i
res = []
unique = 0
for index, row in df2.iterrows():
    temp = [0]*13
    temp[dict1[row['assigned_cluster']]]=1
    res.append(temp)
res = np.array(res)
pkl.dump(res, open("B3Sim/ind.B3Sim1.ally", "wb"))

# Load the cell-cell graph
file1 = open('B3Sim/B3Sim_graph_1.txt', 'r')
Lines = file1.read().splitlines()

# Transform cell-cell graph into .graph file
res = {}
for i in range(dfv.shape[1]):
    res[i]=[]
for i in Lines:
    s = i.split(' ')
    res[int(s[0])].append(int(s[1]))

import pickle as pkl
pkl.dump(res, open("B3Sim/ind.B3Sim1.graph", "wb"))