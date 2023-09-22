## REPLICATE #2 IN SETTY ET AL. DATASET ##

### LIBRARIES ###
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(tidyverse)
library(Azimuth)

### SUBSAMPLING OF REP2 ###
# load processed data from Priyansh
setty_rep2 <- LoadH5Seurat("rep2_Processed_sob.h5seurat")

# subsample processed data
cell.idents <- setty_rep2$orig.ident %>% 
  names %>% enframe(name = NULL, value = "names")

set.seed(4693)
sub.idents <- cell.idents %>% 
  slice_sample(n = 2000)

setty_rep2_sub <- setty_rep2[,sub.idents$names]

# check data structure
DimPlot(setty_rep2_sub, reduction = "pca")
DimPlot(setty_rep2, reduction = "pca")
  
### COMBINING WITH REP1 ###

# load data 
setty_rep1_sub <- LoadH5Seurat("../../HandsOn3_SingleCell_Data_QC/processed_data/setty_seurat_clustered.h5seurat")

# merge with rep1
setty_rep1_sub$orig.ident <- "rep1"

# integrate
setty_list <- list(setty_rep1_sub, 
                   setty_rep2_sub)

features <- SelectIntegrationFeatures(object.list = setty_list)
setty_anchors <- FindIntegrationAnchors(object.list = setty_list,
                                        anchor.features = features)
setty_combined <- IntegrateData(anchorset = setty_anchors)

    # check integration result
    setty_combined <- ScaleData(setty_combined)
    setty_combined <- FindVariableFeatures(setty_combined)
    setty_combined <- RunPCA(setty_combined)
    DimPlot(setty_combined, group.by = "orig.ident")

# clustering
ElbowPlot(setty_combined)

setty_combined <- FindNeighbors(setty_combined, dims = 1:10)
setty_combined <- FindClusters(setty_combined, resolution = 0.5)

  # plot clusters
  DimPlot(setty_combined, 
          reduction = "pca",
          group.by = "seurat_clusters",
          label = TRUE)
  
# assign cell labels
library(bonemarrowref.SeuratData)

setty_combined <- RunAzimuth(query = setty_combined,
                          reference = "bonemarrowref")

  # plot clusters with labels
  setty_combined <- RunTSNE(setty_combined)
  DimPlot(setty_combined, group.by = "predicted.celltype.l2",
          reduction = "tsne", label = TRUE)
  DimPlot(setty_combined, group.by = "orig.ident", 
          reduction = "tsne")
  
# add replicate (simulated) and individual information
setty_combined$individual <- ifelse(setty_combined$orig.ident == "rep1",
                                    yes = "ind1", no = "ind2")
setty_combined$replicate <- c(rep("I1_rep1", 1000),
                              rep("I1_rep2", 1983),
                              rep("I2_rep1", 1000),
                              rep("I2_rep2", 1000))

# save SeuratObject
SaveH5Seurat(setty_combined, 
             filename = "setty_combined_rep1-rep2.h5seurat")


