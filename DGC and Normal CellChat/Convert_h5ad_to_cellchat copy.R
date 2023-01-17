#install.packages("reticulate")
library(reticulate)
library(Seurat)
library(dplyr)
install.packages('NMF')
devtools::install_github("jokergoo/circlize")
devtools::install_github("jokergoo/ComplexHeatmap")

devtools::install_github("sqjin/CellChat")
library(CellChat)
use_condaenv('PY38')

#reticulate::py_discover_config()
#install.packages('anndata')
#anndata::install_anndata()
#library(anndata)

setwd('/Users/gzou/OneDrive - Inside MD Anderson/GSE150290_GCandNT/for Cellchat all pts and 29 normal/final cellchat/')

##############################################################################################################
##################  This is successful only if follow in order ###############################################
##############################################################################################################



### 1. read adata
ad <- import("anndata", convert = FALSE)
########################################################### 2nd try with smaller dataset

ad_object_1 <- ad$read_h5ad("Normal_processed.h5ad")

# access normalized data matrix
data.input <- t(py_to_r(ad_object_1$X))
rownames(data.input) <- rownames(py_to_r(ad_object_1$var))
colnames(data.input) <- rownames(py_to_r(ad_object_1$obs))
# access meta data
meta.data <- py_to_r(ad_object_1$obs)
meta <- meta.data

library(Matrix)
str(data.input)
dgc.Matrix <- as(data.input, "CsparseMatrix")
str(dgc.Matrix)

cellchat <- createCellChat(object = dgc.Matrix, meta = meta, group.by = "leiden_1")
cellchat_1 <- createCellChat(object = dgc.Matrix, meta = meta, group.by = "celltype")

saveRDS(cellchat, 'Normal_cellchat_leiden_1.rds')
saveRDS(cellchat_1, 'Normal_cellchat_celltype.rds')

##############################################################################################
ad_object_2 <- ad$read_h5ad("Type_I_processed.h5ad")

# access normalized data matrix
data.input <- t(py_to_r(ad_object_2$X))
rownames(data.input) <- rownames(py_to_r(ad_object_2$var))
colnames(data.input) <- rownames(py_to_r(ad_object_2$obs))
# access meta data
meta.data <- py_to_r(ad_object_2$obs)
meta <- meta.data

library(Matrix)
str(data.input)
dgc.Matrix <- as(data.input, "CsparseMatrix")
str(dgc.Matrix)

cellchat <- createCellChat(object = dgc.Matrix, meta = meta, group.by = "leiden_1")
cellchat_1 <- createCellChat(object = dgc.Matrix, meta = meta, group.by = "celltype")

saveRDS(cellchat, 'TypeI_cellchat_leiden_1.rds')
saveRDS(cellchat_1, 'TypeI_cellchat_celltype.rds')

##############################################################################################
ad_object_3 <- ad$read_h5ad("Type_II_processed.h5ad")

# access normalized data matrix
data.input <- t(py_to_r(ad_object_3$X))
rownames(data.input) <- rownames(py_to_r(ad_object_3$var))
colnames(data.input) <- rownames(py_to_r(ad_object_3$obs))
# access meta data
meta.data <- py_to_r(ad_object_3$obs)
meta <- meta.data

library(Matrix)
str(data.input)
dgc.Matrix <- as(data.input, "CsparseMatrix")
str(dgc.Matrix)

cellchat <- createCellChat(object = dgc.Matrix, meta = meta, group.by = "leiden_1")
cellchat_1 <- createCellChat(object = dgc.Matrix, meta = meta, group.by = "celltype")

saveRDS(cellchat, 'TypeII_cellchat_leiden_1.rds')
saveRDS(cellchat_1, 'TypeII_cellchat_celltype.rds')
##############################################################################################
ad_object_4 <- ad$read_h5ad("Type_III_processed.h5ad")

# access normalized data matrix
data.input <- t(py_to_r(ad_object_4$X))
rownames(data.input) <- rownames(py_to_r(ad_object_4$var))
colnames(data.input) <- rownames(py_to_r(ad_object_4$obs))
# access meta data
meta.data <- py_to_r(ad_object_4$obs)
meta <- meta.data

library(Matrix)
str(data.input)
dgc.Matrix <- as(data.input, "CsparseMatrix")
str(dgc.Matrix)

cellchat <- createCellChat(object = dgc.Matrix, meta = meta, group.by = "leiden_1")
cellchat_1 <- createCellChat(object = dgc.Matrix, meta = meta, group.by = "celltype")

saveRDS(cellchat, 'TypeIII_cellchat_leiden_1.rds')
saveRDS(cellchat_1, 'TypeIII_cellchat_celltype.rds')