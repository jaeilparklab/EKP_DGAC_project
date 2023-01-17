
library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)
library(ade4)
library(CellChat)

setwd('/Users/gzou/OneDrive - Inside MD Anderson/GSE150290_GCandNT/for Cellchat all pts and 29 normal/final cellchat')


###################################################################

Normal_cellchat.leiden_1 <- readRDS('Normal_cellchat_leiden_1.rds')

levels(Normal_cellchat.leiden_1@idents) # show factor levels of the cell labels

groupSize <- as.numeric(table(Normal_cellchat.leiden_1@idents)) # number of cells in each cell group
#Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling

# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# set the used database in the object
Normal_cellchat.leiden_1@DB <- CellChatDB.use

#Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
Normal_cellchat.leiden_1 <- subsetData(Normal_cellchat.leiden_1) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 8) # do parallel
#> Warning: Strategy 'multiprocess' is deprecated in future (>= 1.20.0). Instead,
#> explicitly specify either 'multisession' or 'multicore'. In the current R
#> session, 'multiprocess' equals 'multisession'.
#> Warning in supportsMulticoreAndRStudio(...): [ONE-TIME WARNING] Forked
#> processing ('multicore') is not supported when running R from RStudio
#> because it is considered unstable. For more details, how to control forked
#> processing or not, and how to silence this warning in future R sessions, see ?
#> parallelly::supportsMulticore
Normal_cellchat.leiden_1 <- identifyOverExpressedGenes(Normal_cellchat.leiden_1)
Normal_cellchat.leiden_1 <- identifyOverExpressedInteractions(Normal_cellchat.leiden_1)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
#Normal_cellchat.leiden_1 <- projectData(Normal_cellchat.leiden_1, PPI.human)
#Part II: Inference of cell-cell communication network
#Compute the communication probability and infer cellular communication network

Normal_cellchat.leiden_1 <- computeCommunProb(Normal_cellchat.leiden_1)

#Normal_cellchat.leiden_1 <- computeCommunProb(
#  Normal_cellchat.leiden_1,
#  type = c("triMean", "truncatedMean", "median"),
#  trim = NULL,
#  LR.use = NULL,
#  raw.use = TRUE,
#  population.size = FALSE,
#  do.fast = TRUE,
#  nboot = 100,
#  seed.use = 1L,
#  Kh = 0.5,
#  n = 1
#)


# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
Normal_cellchat.leiden_1 <- filterCommunication(Normal_cellchat.leiden_1, min.cells = 10)

#Extract the inferred cellular communication network as a data frame
#table(Normal_cellchat.leiden_1@meta$group)
#df.net <- subsetCommunication(Normal_cellchat.leiden_1, sources.use = c(1,13,16,17,19), targets.use = c(2,3,4,5,6,7,8,9,10,11,12,14,15,18,20,21,22)) #gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.

#df.net <- subsetCommunication(Normal_cellchat.leiden_1, signaling = c("CCL")) #gives the inferred cell-cell communications mediated by signaling WNT and TGFb.
#Infer the cell-cell communication at a signaling pathway level
Normal_cellchat.leiden_1 <- computeCommunProbPathway(Normal_cellchat.leiden_1)
#Calculate the aggregated cell-cell communication network
Normal_cellchat.leiden_1 <- aggregateNet(Normal_cellchat.leiden_1)

groupSize <- as.numeric(table(Normal_cellchat.leiden_1@idents))
par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(Normal_cellchat.leiden_1@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(Normal_cellchat.leiden_1@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


mat <- Normal_cellchat.leiden_1@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

#Visualize each signaling pathway using Hierarchy plot, Circle plot or Chord diagram
Normal_cellchat.leiden_1@netP$pathways

pathways.show <- c("MHC-I") #"CXCL","PTN","MHC-I","NOTCH"
# Hierarchy plot


#[1] "B cell"          "Epithelial cell" "Mast cell"       "Myeloid cell"    "T cell" 
vertex.receiver = c(1,2,3,4,7,8,10,12,13,14,18) # a numeric vector. 
netVisual_aggregate(Normal_cellchat.leiden_1, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(Normal_cellchat.leiden_1, signaling = pathways.show, layout = "circle")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(Normal_cellchat.leiden_1, signaling = pathways.show, layout = "chord")
#[1] "B cell"          "Epithelial cell" "Fibroblast"      "Mast cell"       "Myeloid cell"    "T cell"

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')


# show all the interactions between epithelial cell and T cell 
netVisual_chord_gene(Normal_cellchat.leiden_1, sources.use = c(1,2,3,4,7,8,10,12,13,14,18), targets.use = c(5,6,9), lab.cex = 0.5,legend.pos.y = 30)

# show all the interactions between epithelial cell and Myeloid 
netVisual_chord_gene(Normal_cellchat.leiden_1, sources.use = c(1,2,3,4,7,8,10,12,13,14,18), targets.use = c(11,19), lab.cex = 0.5,legend.pos.y = 30)

# show all the interactions between epithelial cell and Bcell 
netVisual_chord_gene(Normal_cellchat.leiden_1, sources.use = c(1,2,3,4,7,8,10,12,13,14,18), targets.use = c(16,17), lab.cex = 0.5,legend.pos.y = 30)

# show all the interactions between epithelial cells and fibroblast
netVisual_chord_gene(Normal_cellchat.leiden_1, sources.use = c(1,2,3,4,7,8,10,12,13,14,18), targets.use = c(15), lab.cex = 0.5,legend.pos.x = 15,legend.pos.y = 10)


# show all the significant interactions (L-R pairs) associated with certain signaling pathways
#epithelial cell and T cell
netVisual_bubble(Normal_cellchat.leiden_1, sources.use = c(1,2,3,4,7,8,10,12,13,14,18), targets.use = c(5,6,9), signaling = c("MHC-I","MHC-II"), remove.isolate = FALSE)
#> Comparing communications on a single object
saveRDS(Normal_cellchat.leiden_1, file = 'Normal_cellchat.leiden_1_after_cellchat.rds')





###################################################################

IP_007_cellchat.celltype <- readRDS('IP_007_cellchat_celltype.rds')
levels(IP_007_cellchat.celltype@idents) # show factor levels of the cell labels

groupSize <- as.numeric(table(IP_007_cellchat.celltype@idents)) # number of cells in each cell group
#Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling

# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# set the used database in the object
IP_007_cellchat.celltype@DB <- CellChatDB.use

#Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
IP_007_cellchat.celltype <- subsetData(IP_007_cellchat.celltype) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 8) # do parallel
#> Warning: Strategy 'multiprocess' is deprecated in future (>= 1.20.0). Instead,
#> explicitly specify either 'multisession' or 'multicore'. In the current R
#> session, 'multiprocess' equals 'multisession'.
#> Warning in supportsMulticoreAndRStudio(...): [ONE-TIME WARNING] Forked
#> processing ('multicore') is not supported when running R from RStudio
#> because it is considered unstable. For more details, how to control forked
#> processing or not, and how to silence this warning in future R sessions, see ?
#> parallelly::supportsMulticore
IP_007_cellchat.celltype <- identifyOverExpressedGenes(IP_007_cellchat.celltype)
IP_007_cellchat.celltype <- identifyOverExpressedInteractions(IP_007_cellchat.celltype)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
#IP_007_cellchat.celltype <- projectData(IP_007_cellchat.celltype, PPI.human)
#Part II: Inference of cell-cell communication network
#Compute the communication probability and infer cellular communication network

IP_007_cellchat.celltype <- computeCommunProb(IP_007_cellchat.celltype)

#IP_007_cellchat.celltype <- computeCommunProb(
#  IP_007_cellchat.celltype,
#  type = c("triMean", "truncatedMean", "median"),
#  trim = NULL,
#  LR.use = NULL,
#  raw.use = TRUE,
#  population.size = FALSE,
#  do.fast = TRUE,
#  nboot = 100,
#  seed.use = 1L,
#  Kh = 0.5,
#  n = 1
#)


# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
IP_007_cellchat.celltype <- filterCommunication(IP_007_cellchat.celltype, min.cells = 10)

#Extract the inferred cellular communication network as a data frame
#table(IP_007_cellchat.celltype@meta$group)

#[1] "Epi_1"      "Tcell_1"    "B cell_1"   "Macro_1"    "Monocyte_1" "Tcell_3"    "Tcell_4"    "Tcell_5"    "Tcell_6"    "Monocyte_2" "Monocyte_3" "Monocyte_4"
#[13] "Epi_2"      "Bcell_2"    "Mast cell"  "Epi_3"      "Epi_4"      "Macro_2"    "Epi_6"      "Monocyte_5" "Bcell_3"    "Tcell_8" 

df.net <- subsetCommunication(IP_007_cellchat.celltype, sources.use = c(1,13,16,17,19), targets.use = c(2,3,4,5,6,7,8,9,10,11,12,14,15,18,20,21,22)) #gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.

df.net <- subsetCommunication(IP_007_cellchat.celltype, signaling = c("CCL")) #gives the inferred cell-cell communications mediated by 
#Infer the cell-cell communication at a signaling pathway level
IP_007_cellchat.celltype <- computeCommunProbPathway(IP_007_cellchat.celltype)
#Calculate the aggregated cell-cell communication network
IP_007_cellchat.celltype <- aggregateNet(IP_007_cellchat.celltype)

groupSize <- as.numeric(table(IP_007_cellchat.celltype@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(IP_007_cellchat.celltype@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(IP_007_cellchat.celltype@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


mat <- IP_007_cellchat.celltype@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

#Visualize each signaling pathway using Hierarchy plot, Circle plot or Chord diagram
IP_007_cellchat.celltype@netP$pathways
#[1] "MIF"      "MHC-I"    "COLLAGEN" "CD99"     "APP"      "CXCL"     "CLEC"     "LAMININ"  "MPZ"      "FN1"      "THBS"     "VISFATIN" "FGF"      "IL1"      "TIGIT"   
#[16] "NECTIN"   "PECAM1"   "CALCR"    "ANGPT"    "MK"       "GRN"      "GAS"      "EGF"      "KIT"      "NOTCH"    "CD46"     "JAM"      "PVR"      "CD96"   
pathways.show <- c("CCL",'CXCL' ) #"CXCL","PTN","MHC-I","NOTCH"
# Hierarchy plot


#[1] "Epi_1"      "Tcell_1"    "B cell_1"   "Macro_1"    "Monocyte_1" "Tcell_3"    "Tcell_4"    "Tcell_5"    "Tcell_6"    "Monocyte_2" "Monocyte_3" "Monocyte_4"
#[13] "Epi_2"      "Bcell_2"    "Mast cell"  "Epi_3"      "Epi_4"      "Macro_2"    "Epi_6"      "Monocyte_5" "Bcell_3"    "Tcell_8" 
vertex.receiver = seq(2,4,5) # a numeric vector. 
netVisual_aggregate(IP_007_cellchat.celltype, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(IP_007_cellchat.celltype, signaling = pathways.show, layout = "circle")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(IP_007_cellchat.celltype, signaling = pathways.show, layout = "chord")

#[1] "Epi_1"      "Tcell_1"    "B cell_1"   "Macro_1"    "Monocyte_1" "Tcell_3"    "Tcell_4"    "Tcell_5"    "Tcell_6"    "Monocyte_2" "Monocyte_3" "Monocyte_4"
#[13] "Epi_2"      "Bcell_2"    "Mast cell"  "Epi_3"      "Epi_4"      "Macro_2"    "Epi_6"      "Monocyte_5" "Bcell_3"    "Tcell_8" 

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')


# show all the interactions between epithelial cell and T cell 
netVisual_chord_gene(IP_007_cellchat.celltype, sources.use = c(1,13,16,17,19), targets.use = c(2,6,7,8,9,22), lab.cex = 0.5,legend.pos.y = 30)

# show all the interactions between epithelial cell and Monocytes 
netVisual_chord_gene(IP_007_cellchat.celltype, sources.use = c(1,13,16,17,19), targets.use = c(5,10,11,12,20), lab.cex = 0.5,legend.pos.y = 30)

# show all the interactions between epithelial cells and [T cell and Monocyte]
netVisual_chord_gene(IP_007_cellchat.celltype, sources.use = c(1,13,16,17,19), targets.use = c(2,6,7,8,9,22,5,10,11,12,20), lab.cex = 0.5,legend.pos.x = 15,legend.pos.y = 10)


# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_chord_gene(IP_007_cellchat.celltype, sources.use = c(1,13,16,17,19), targets.use = c(2,6,7,8,9,22,5,10,11,12,20),  signaling = c("CCL","CXCL"), lab.cex = 0.5,legend.pos.x = 15,legend.pos.y = 10)
netVisual_bubble(P009, sources.use = c(1,13,16,17,19), targets.use = c(2,6,7,8,9,22,5,10,11,12,20), signaling = c("CCL","CXCL"), remove.isolate = FALSE)
#> Comparing communications on a single object

saveRDS(IP_007_cellchat.celltype, file = 'IP_007_cellchat.celltype_after_cellchat.rds')
