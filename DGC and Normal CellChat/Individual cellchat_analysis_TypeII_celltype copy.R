
library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)
library(ade4)
library(CellChat)

setwd('/Users/gzou/OneDrive - Inside MD Anderson/GSE150290_GCandNT/for Cellchat all pts and 29 Normal/final cellchat')


###################################################################

TypeII_cellchat.celltype <- readRDS('TypeII_cellchat_celltype.rds')

levels(TypeII_cellchat.celltype@idents) # show factor levels of the cell labels

groupSize <- as.numeric(table(TypeII_cellchat.celltype@idents)) # number of cells in each cell group
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
TypeII_cellchat.celltype@DB <- CellChatDB.use

#Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
TypeII_cellchat.celltype <- subsetData(TypeII_cellchat.celltype) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 8) # do parallel
#> Warning: Strategy 'multiprocess' is deprecated in future (>= 1.20.0). Instead,
#> explicitly specify either 'multisession' or 'multicore'. In the current R
#> session, 'multiprocess' equals 'multisession'.
#> Warning in supportsMulticoreAndRStudio(...): [ONE-TIME WARNING] Forked
#> processing ('multicore') is not supported when running R from RStudio
#> because it is considered unstable. For more details, how to control forked
#> processing or not, and how to silence this warning in future R sessions, see ?
#> parallelly::supportsMulticore
TypeII_cellchat.celltype <- identifyOverExpressedGenes(TypeII_cellchat.celltype)
TypeII_cellchat.celltype <- identifyOverExpressedInteractions(TypeII_cellchat.celltype)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
#TypeII_cellchat.celltype <- projectData(TypeII_cellchat.celltype, PPI.human)
#Part II: Inference of cell-cell communication network
#Compute the communication probability and infer cellular communication network

TypeII_cellchat.celltype <- computeCommunProb(TypeII_cellchat.celltype)

#TypeII_cellchat.celltype <- computeCommunProb(
#  TypeII_cellchat.celltype,
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
TypeII_cellchat.celltype <- filterCommunication(TypeII_cellchat.celltype, min.cells = 10)

#Extract the inferred cellular communication network as a data frame
#table(TypeII_cellchat.celltype@meta$group)
#df.net <- subsetCommunication(TypeII_cellchat.celltype, sources.use = c(1,13,16,17,19), targets.use = c(2,3,4,5,6,7,8,9,10,11,12,14,15,18,20,21,22)) #gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.

#df.net <- subsetCommunication(TypeII_cellchat.celltype, signaling = c("CCL")) #gives the inferred cell-cell communications mediated by signaling WNT and TGFb.
#Infer the cell-cell communication at a signaling pathway level
TypeII_cellchat.celltype <- computeCommunProbPathway(TypeII_cellchat.celltype)
#Calculate the aggregated cell-cell communication network
TypeII_cellchat.celltype <- aggregateNet(TypeII_cellchat.celltype)

groupSize <- as.numeric(table(TypeII_cellchat.celltype@idents))
par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(TypeII_cellchat.celltype@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(TypeII_cellchat.celltype@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


mat <- TypeII_cellchat.celltype@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

#Visualize each signaling pathway using Hierarchy plot, Circle plot or Chord diagram
TypeII_cellchat.celltype@netP$pathways

pathways.show <- c("MHC-I") #"CXCL","PTN","MHC-I","NOTCH"
# Hierarchy plot


#[1] "B cell"          "Epithelial cell" "Mast cell"       "Myeloid cell"    "T cell" 
vertex.receiver = c(1,2,3,4,7,8,10,12,13,14,18) # a numeric vector. 
netVisual_aggregate(TypeII_cellchat.celltype, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(TypeII_cellchat.celltype, signaling = pathways.show, layout = "circle")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(TypeII_cellchat.celltype, signaling = pathways.show, layout = "chord")
#[1] "B cell"          "Epithelial cell" "Fibroblast"      "Mast cell"       "Myeloid cell"    "T cell"

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')


# show all the interactions between epithelial cell and T cell 
netVisual_chord_gene(TypeII_cellchat.celltype, sources.use = c(1,2,3,4,7,8,10,12,13,14,18), targets.use = c(5,6,9), lab.cex = 0.5,legend.pos.y = 30)

# show all the interactions between epithelial cell and Myeloid 
netVisual_chord_gene(TypeII_cellchat.celltype, sources.use = c(1,2,3,4,7,8,10,12,13,14,18), targets.use = c(11,19), lab.cex = 0.5,legend.pos.y = 30)

# show all the interactions between epithelial cell and Bcell 
netVisual_chord_gene(TypeII_cellchat.celltype, sources.use = c(1,2,3,4,7,8,10,12,13,14,18), targets.use = c(16,17), lab.cex = 0.5,legend.pos.y = 30)

# show all the interactions between epithelial cells and fibroblast
netVisual_chord_gene(TypeII_cellchat.celltype, sources.use = c(1,2,3,4,7,8,10,12,13,14,18), targets.use = c(15), lab.cex = 0.5,legend.pos.x = 15,legend.pos.y = 10)


# show all the significant interactions (L-R pairs) associated with certain signaling pathways
#epithelial cell and T cell
netVisual_bubble(TypeII_cellchat.celltype, sources.use = c(1,2,3,4,7,8,10,12,13,14,18), targets.use = c(5,6,9), signaling = c("MHC-I","MHC-II"), remove.isolate = FALSE)
#> Comparing communications on a single object
saveRDS(TypeII_cellchat.celltype, file = 'TypeII_cellchat.celltype_after_cellchat.rds')
