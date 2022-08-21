rm(list = ls())
if(!is.null(dev.list())) dev.off()
cat("\f")

# The code is based on the Seurat algorithm (Stuart et al., 2019)
# The source code can be found at https://github.com/satijalab/seurat

library(dplyr)
library(Seurat)
library(textshape)


# Exemplery clustering codes for 1694-cell Simulated datasets
for (i in 1:10){
  # Load one of ten simulated datasets
  D1200.sci.log<-read.csv(paste0("/Users/lbyjo/Desktop/Pingzhao Hu/Simulation/DSim/DSim_",i,".csv"), row.names = 1)
  
  # Pre-processing
  D <- CreateSeuratObject(counts = D1200.sci.log)
  D <- FindVariableFeatures(D, nfeatures = 1200)
  all.genes.D <- rownames(D)
  D <- ScaleData(D, features = all.genes.D)
  
  # Clustering
  D <- RunPCA(D, features = VariableFeatures(object = D))
  D <- FindNeighbors(D, dims = 1:10)
  D <- FindClusters(D, resolution = 0.5)
  saveRDS(Idents(D), paste0('/Users/lbyjo/Desktop/Pingzhao Hu/Simulation/RDS/DSim_',i,'.RDS'))
}



