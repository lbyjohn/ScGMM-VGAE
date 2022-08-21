rm(list = ls())
if(!is.null(dev.list())) dev.off()
cat("\f")

# Read raw data
Baron3<-read.csv("/Users/lbyjo/Desktop/Pingzhao Hu/Rawdata/Baron3.csv")
Baron4<-read.csv("/Users/lbyjo/Desktop/Pingzhao Hu/Rawdata/Baron4.csv")
Darmanis<-read.csv("/Users/lbyjo/Desktop/Pingzhao Hu/Rawdata/Darmanis.csv", row.names = 1)


# Transpose datasets (rows:genes, columns:cell)
Baron3.t<-Baron3[,4:20128]
Baron3.t<-data.frame(t(Baron3.t))
colnames(Baron3.t)<-Baron3$barcode

Baron4.t<-Baron4[,4:20128]
Baron4.t<-data.frame(t(Baron4.t))
colnames(Baron4.t)<-Baron4$barcode


# Log-Normalize the data
library(Seurat)
D.norm.log <- NormalizeData(Darmanis, normalization.method = "LogNormalize", scale.factor = 10000)
D.norm.log<-data.frame(D.norm.log)
write.csv(D.norm.log,"/Users/lbyjo/Desktop/Pingzhao Hu/Impute/norm_data/D_norm_log.csv", row.names = TRUE)

B3.norm.log <- NormalizeData(Baron3.t, normalization.method = "LogNormalize", scale.factor = 10000)
B3.norm.log<-data.frame(B3.norm.log)
write.csv(B3.norm.log,"/Users/lbyjo/Desktop/Pingzhao Hu/Impute/norm_data/B3_norm_log.csv", row.names = TRUE)

B4.norm.log <- NormalizeData(Baron4.t, normalization.method = "LogNormalize", scale.factor = 10000)
B4.norm.log<-data.frame(B4.norm.log)
write.csv(B4.norm.log,"/Users/lbyjo/Desktop/Pingzhao Hu/Impute/norm_data/B4_norm_log.csv", row.names = TRUE)


# Load the imputation methods
library(scImpute)

# Imputation
# ScImpute
scimpute(
  count_path="/Users/lbyjo/Desktop/Pingzhao Hu/Impute/norm_data/D_norm_log.csv",
  Kcluster=8, out_dir ="/Users/lbyjo/Desktop/Pingzhao Hu/Impute/Scimpute_D",ncores=1)
D.sci<-read.csv("/Users/lbyjo/Desktop/Pingzhao Hu/Impute/Scimpute_D/Scimpute_Dscimpute_count.csv",row.names=1)

scimpute(
  count_path="/Users/lbyjo/Desktop/Pingzhao Hu/Impute/norm_data/B3_norm_log.csv",
  Kcluster=14, out_dir ="/Users/lbyjo/Desktop/Pingzhao Hu/Impute/Scimpute_B3",ncores=1)
B3.sci<-read.csv("/Users/lbyjo/Desktop/Pingzhao Hu/Impute/Scimpute_B3/Scimpute_B3scimpute_count.csv",row.names=1)

scimpute(
  count_path="/Users/lbyjo/Desktop/Pingzhao Hu/Impute/norm_data/B4_norm_log.csv",
  Kcluster=14, out_dir ="/Users/lbyjo/Desktop/Pingzhao Hu/Impute/Scimpute_B4",ncores=1)
B4.sci<-read.csv("/Users/lbyjo/Desktop/Pingzhao Hu/Impute/Scimpute_B4/Scimpute_B4scimpute_count.csv",row.names=1)


# Select the top 1200 genes
#ScImpute
D.sci_seurat <- CreateSeuratObject(counts = D.sci)
D.sci.selected <- FindVariableFeatures(D.sci_seurat, nfeatures = 1200)
D.sci.topfeat = D.sci.selected@assays[["RNA"]]@var.features
D.sci.selected.data <- as.data.frame(GetAssayData(object = D.sci.selected))[D.sci.topfeat, ]
write.csv(D.sci.selected.data,"/Users/lbyjo/Desktop/Pingzhao Hu/Impute/top1200/D_sci_Top1200.csv", row.names = TRUE)

B3.sci_seurat <- CreateSeuratObject(counts = B3.sci)
B3.sci.selected <- FindVariableFeatures(B3.sci_seurat, nfeatures = 1200)
B3.sci.topfeat = B3.sci.selected@assays[["RNA"]]@var.features
B3.sci.selected.data <- as.data.frame(GetAssayData(object = B3.sci.selected))[B3.sci.topfeat, ]
write.csv(B3.sci.selected.data,"/Users/lbyjo/Desktop/Pingzhao Hu/Impute/top1200/B3_sci_Top1200.csv", row.names = TRUE)

B4.sci_seurat <- CreateSeuratObject(counts = B4.sci)
B4.sci.selected <- FindVariableFeatures(B4.sci_seurat, nfeatures = 1200)
B4.sci.topfeat = B4.sci.selected@assays[["RNA"]]@var.features
B4.sci.selected.data <- as.data.frame(GetAssayData(object = B4.sci.selected))[B4.sci.topfeat, ]
write.csv(B4.sci.selected.data,"/Users/lbyjo/Desktop/Pingzhao Hu/Impute/top1200/B4_sci_Top1200.csv", row.names = TRUE)


# Generate cell-cell graph
library(SingleCellExperiment)
library(scRNAseq)
library(scran)
library(igraph)
library(BiocSingular)
# Read the filtered data
bsparam()
options(BiocSingularParam.default=ExactParam())
bsparam()

D.sci.graph <- buildKNNGraph(D.sci.selected.data, k=5, d=10)
write_graph(D.sci.graph, "/Users/lbyjo/Desktop/Pingzhao Hu/Impute/cellgraph/D_sci_graph.txt", format = c("edgelist"))
B3.sci.graph <- buildKNNGraph(B3.sci.selected.data, k=5, d=10)
write_graph(B3.sci.graph, "/Users/lbyjo/Desktop/Pingzhao Hu/Impute/cellgraph/B3_sci_graph.txt", format = c("edgelist"))
B4.sci.graph <- buildKNNGraph(B4.sci.selected.data, k=5, d=10)
write_graph(B4.sci.graph, "/Users/lbyjo/Desktop/Pingzhao Hu/Impute/cellgraph/B4_sci_graph.txt", format = c("edgelist"))





