rm(list = ls())
if(!is.null(dev.list())) dev.off()
cat("\f")


########### Simulation using Zheng preset
# Load raw dataset
library(Seurat)
raw=ReadMtx(mtx='/Users/lbyjo/Desktop/Pingzhao Hu/Rawdata/filtered_matrices_mex/hg19/matrix.mtx', 
            cells = '/Users/lbyjo/Desktop/Pingzhao Hu/Rawdata/filtered_matrices_mex/hg19/barcodes.tsv',
            features = '/Users/lbyjo/Desktop/Pingzhao Hu/Rawdata/filtered_matrices_mex/hg19/genes.tsv')
raw=data.frame(raw)

# The raw data contains 32738 genes, while the template for simulation contains 19536 genes
# So select those 19536 genes from the raw data.
library(SPARSim)
data(Zheng_param_preset)
selectgenes<-attributes(Zheng_param_preset$Zheng_C1$intensity)$names
raw_filter<-raw[selectgenes,]
raw_filter[is.na(raw_filter)]<-0

# Filter 1200 genes
Z_seurat <- CreateSeuratObject(counts = raw_filter)
Z_lognorm <- NormalizeData(Z_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
Z.selected <- FindVariableFeatures(Z_lognorm, nfeatures = 1200)
top1200 = Z.selected@assays[["RNA"]]@var.features

# Load the parameter presets for each cluster
Zheng_C1<-Zheng_param_preset$Zheng_C1
Zheng_C2<-Zheng_param_preset$Zheng_C2
Zheng_C3<-Zheng_param_preset$Zheng_C3
Zheng_C4<-Zheng_param_preset$Zheng_C4

# Load the cell graph generation package
library(SingleCellExperiment)
library(scRNAseq)
library(scran)
library(igraph)
library(BiocSingular)
bsparam()
options(BiocSingularParam.default=ExactParam())
bsparam()


#### Simulate same-size dataset (3388 cells)
# Change the parameters to only simulate the top1200 genes
cond_A_param <- SPARSim_create_simulation_parameter(
  intensity = Zheng_C1$intensity[top1200], 
  variability = Zheng_C1$variability[top1200], 
  library_size = Zheng_C1$lib_size, 
  condition_name = "cond_A")

cond_B_param <- SPARSim_create_simulation_parameter(
  intensity = Zheng_C2$intensity[top1200], 
  variability = Zheng_C2$variability[top1200], 
  library_size = Zheng_C2$lib_size, 
  condition_name = "cond_B")

cond_C_param <- SPARSim_create_simulation_parameter(
  intensity = Zheng_C3$intensity[top1200], 
  variability = Zheng_C3$variability[top1200], 
  library_size = Zheng_C3$lib_size, 
  condition_name = "cond_C")

cond_D_param <- SPARSim_create_simulation_parameter(
  intensity = Zheng_C4$intensity[top1200], 
  variability = Zheng_C4$variability[top1200], 
  library_size = Zheng_C4$lib_size, 
  condition_name = "cond_D")

SPARSim_sim_param <- list(cond_A = cond_A_param, cond_B = cond_B_param, cond_C = cond_C_param, cond_D = cond_D_param)

for (i in 1:10){
  set.seed(i)
  # Simulate dataset
  sim_result <- SPARSim_simulation(dataset_parameter = SPARSim_sim_param)
  Zheng_same<-data.frame(sim_result$count_matrix)
  write.csv(Zheng_same,paste0("/Users/lbyjo/Desktop/Pingzhao Hu/Simulation/test/Z_same_",i,".csv"), row.names = TRUE)
  # Generate cell graph
  Z.graph <- buildKNNGraph(Zheng_same, k=5, d=10)
  write_graph(Z.graph, paste0("/Users/lbyjo/Desktop/Pingzhao Hu/Simulation/test/Z_same_graph_",i,".txt"), format = c("edgelist"))
}

# Create Labels
cluster_same<-c(1440,1718,184,46)
Z_same_label<-Zheng_same[1,]
for (i in 1:4){
  Z_same_label[1,(sum(cluster_same[0:(i-1)])+1):sum(cluster_same[0:i])]<-i
}
row.names(Z_same_label)<-'assigned_cluster' 
Z_same_label<-data.frame(t(Z_same_label))
write.csv(Z_same_label,"/Users/lbyjo/Desktop/Pingzhao Hu/Simulation/test/Z_same_label.csv", row.names = TRUE)



##### Simulate half-size dataset (1694 cells)
# Change the parameters to only simulate the top1200 genes
set.seed(123)
cond_A_param <- SPARSim_create_simulation_parameter(
  intensity = Zheng_C1$intensity[top1200], 
  variability = Zheng_C1$variability[top1200], 
  library_size = sample(Zheng_C1$lib_size, size =  720, replace = FALSE), 
  condition_name = "cond_A")

cond_B_param <- SPARSim_create_simulation_parameter(
  intensity = Zheng_C2$intensity[top1200], 
  variability = Zheng_C2$variability[top1200], 
  library_size = sample(Zheng_C2$lib_size, size =  859, replace = FALSE), 
  condition_name = "cond_B")

cond_C_param <- SPARSim_create_simulation_parameter(
  intensity = Zheng_C3$intensity[top1200], 
  variability = Zheng_C3$variability[top1200], 
  library_size = sample(Zheng_C3$lib_size, size =  92, replace = FALSE), 
  condition_name = "cond_C")

cond_D_param <- SPARSim_create_simulation_parameter(
  intensity = Zheng_C4$intensity[top1200], 
  variability = Zheng_C4$variability[top1200], 
  library_size = sample(Zheng_C4$lib_size, size =  23, replace = FALSE), 
  condition_name = "cond_D")

SPARSim_sim_param <- list(cond_A = cond_A_param, cond_B = cond_B_param, cond_C = cond_C_param, cond_D = cond_D_param)

for (i in 1:10){
  set.seed(i)
  # Simulate dataset
  sim_result <- SPARSim_simulation(dataset_parameter = SPARSim_sim_param)
  Zheng_half<-data.frame(sim_result$count_matrix)
  write.csv(Zheng_half,paste0("/Users/lbyjo/Desktop/Pingzhao Hu/Simulation/test/Z_half_",i,".csv"), row.names = TRUE)
  # Generate cell graph
  Z.graph <- buildKNNGraph(Zheng_half, k=5, d=10)
  write_graph(Z.graph, paste0("/Users/lbyjo/Desktop/Pingzhao Hu/Simulation/test/Z_half_graph_",i,".txt"), format = c("edgelist"))
}

# Create Labels
cluster_half<-c(720, 859,92,23)
Z_half_label<-Zheng_half[1,]
for (i in 1:4){
  Z_half_label[1,(sum(cluster_half[0:(i-1)])+1):sum(cluster_half[0:i])]<-i
}
row.names(Z_half_label)<-'assigned_cluster' 
Z_half_label<-data.frame(t(Z_half_label))
write.csv(Z_half_label,"/Users/lbyjo/Desktop/Pingzhao Hu/Simulation/test/Z_half_label.csv", row.names = TRUE)



##### Simulate mid-size dataset (6776 cells)
# Change the parameters to only simulate the top1200 genes
set.seed(123)
cond_A_param <- SPARSim_create_simulation_parameter(
  intensity = Zheng_C1$intensity[top1200], 
  variability = Zheng_C1$variability[top1200], 
  library_size = sample(Zheng_C1$lib_size, size =  2880, replace = TRUE), 
  condition_name = "cond_A")

cond_B_param <- SPARSim_create_simulation_parameter(
  intensity = Zheng_C2$intensity[top1200], 
  variability = Zheng_C2$variability[top1200], 
  library_size = sample(Zheng_C2$lib_size, size =  3436, replace = TRUE), 
  condition_name = "cond_B")

cond_C_param <- SPARSim_create_simulation_parameter(
  intensity = Zheng_C3$intensity[top1200], 
  variability = Zheng_C3$variability[top1200], 
  library_size = sample(Zheng_C3$lib_size, size =  368, replace = TRUE), 
  condition_name = "cond_C")

cond_D_param <- SPARSim_create_simulation_parameter(
  intensity = Zheng_C4$intensity[top1200], 
  variability = Zheng_C4$variability[top1200], 
  library_size = sample(Zheng_C4$lib_size, size =  92, replace = TRUE), 
  condition_name = "cond_D")

SPARSim_sim_param <- list(cond_A = cond_A_param, cond_B = cond_B_param, cond_C = cond_C_param, cond_D = cond_D_param)

for (i in 1:10){
  set.seed(i)
  # Simulate dataset
  sim_result <- SPARSim_simulation(dataset_parameter = SPARSim_sim_param)
  Zheng_mid<-data.frame(sim_result$count_matrix)
  write.csv(Zheng_mid,paste0("/Users/lbyjo/Desktop/Pingzhao Hu/Simulation/test/Z_mid_",i,".csv"), row.names = TRUE)
  # Generate cell graph
  Z.graph <- buildKNNGraph(Zheng_mid, k=5, d=10)
  write_graph(Z.graph, paste0("/Users/lbyjo/Desktop/Pingzhao Hu/Simulation/test/Z_mid_graph_",i,".txt"), format = c("edgelist"))
}

# Create Labels
cluster_mid<-c(2880, 3436,368,92)
Z_mid_label<-Zheng_mid[1,]
for (i in 1:4){
  Z_mid_label[1,(sum(cluster_mid[0:(i-1)])+1):sum(cluster_mid[0:i])]<-i
}
row.names(Z_mid_label)<-'assigned_cluster' 
Z_mid_label<-data.frame(t(Z_mid_label))
write.csv(Z_mid_label,"/Users/lbyjo/Desktop/Pingzhao Hu/Simulation/test/Z_mid_label.csv", row.names = TRUE)






