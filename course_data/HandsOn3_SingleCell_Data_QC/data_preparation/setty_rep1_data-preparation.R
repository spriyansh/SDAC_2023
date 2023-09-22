## REPLICATE #1 IN SETTY ET AL. DATASET ##

### LIBRARIES ###
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(tidyverse)

# load raw data
setty_rep1_matrix <- ReadMtx(mtx = "../../HandsOn2_OmicsBox_Count_Tables_Filtered/rep1_p1_downSampled/Rep1_P1_DownSampled_OB_Counts.mtx",
                             cells = "../../HandsOn2_OmicsBox_Count_Tables_Filtered/rep1_p1_downSampled/Rep1_P1_DownSampled_OB_Counts_barcodes.tsv.gz",
                             features = "../../HandsOn2_OmicsBox_Count_Tables_Filtered/rep1_p1_downSampled/Rep1_P1_DownSampled_OB_Counts_features.tsv.gz")
setty_rep1 <- CreateSeuratObject(setty_rep1_matrix,
                                 project = "healthyBM",
                                 min.cells = 1)

# subsample raw data
cell.idents <- setty_rep1$orig.ident %>% 
  names %>% enframe(name = NULL, value = "names")

set.seed(4693)
sub.idents <- cell.idents %>% 
  slice_sample(n = 3000)

setty_rep1_sub <- setty_rep1[,sub.idents$names]

# check data structure

  # subsampled
  setty_rep1_sub <- FindVariableFeatures(setty_rep1_sub, 
                                         selection.method = "vst",
                                         nfeatures = 2000)
  setty_rep1_sub <- ScaleData(setty_rep1_sub)
  setty_rep1_sub <- RunPCA(setty_rep1_sub)
  DimPlot(setty_rep1_sub)
  
  # full
  setty_rep1 <- FindVariableFeatures(setty_rep1, 
                                         selection.method = "vst",
                                         nfeatures = 2000)
  setty_rep1 <- ScaleData(setty_rep1)
  setty_rep1 <- RunPCA(setty_rep1)
  DimPlot(setty_rep1)
  
# export seurat file

  # reset
  setty_rep1 <- CreateSeuratObject(setty_rep1_matrix,
                                   project = "healthyBM",
                                   min.cells = 1)
  setty_rep1_sub <- setty_rep1[,sub.idents$names]
  
  # save
  SaveH5Seurat(setty_rep1_sub, filename = "setty_seurat.h5seurat")
  
  
