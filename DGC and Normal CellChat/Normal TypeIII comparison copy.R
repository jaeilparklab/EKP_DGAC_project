##################################################################################################################
#####merge Normal.celltype and TypeIII.celltype---merge Normal.celltype and TypeIII.celltype######################################
#####merge Normal.celltype and TypeIII.celltype---merge Normal.celltype and TypeIII.celltype######################################
#####merge Normal.celltype and TypeIII.celltype---merge Normal.celltype and TypeIII.celltype######################################
##WT KO comparison from here################################################################################################################
library(CellChat)
library(ggplot2)                  
library(patchwork)
library(igraph)
setwd('/Users/gzou/OneDrive - Inside MD Anderson/GSE150290_GCandNT/for Cellchat all pts and 29 normal/Cellchat/')

Normal.celltype <- readRDS("Normal_cellchat.celltype_after_cellchat.rds")
Normal.celltype
TypeIII.celltype <- readRDS("TypeIII_cellchat.celltype_after_cellchat.rds")
TypeIII.celltype

levels(Normal.celltype@idents)
levels(TypeIII.celltype@idents)
#Lift up CellChat object and merge together
#Since there are additional two populations (i.e., dermal DC and pericytes) specific to E14.5 
#compared to E13.5, we lift up cellchat.E13 by lifting up the cell groups to the same cell 
#labels as E14.5. liftCellChat will only update the slot related to the cell-cell 
#communication network, including slots object@net, object@netP and object@idents.

# Define the cell labels to lift up, if the idents in one group can cover the other group
#group.new = levels(TypeIII.celltype@idents)
#Normal.celltype <- liftCellChat(Normal.celltype, group.new)
#levels(Normal.celltype@idents)

#> The CellChat object will be lifted up using the cell labels FIB-A, FIB-B, FIB-P, DC, Pericyte, MYL, Immune, ENDO, Muscle, MELA, Basal-P, Basal, Spinious
#> Update slots object@net, object@netP, object@idents in a single dataset...
#object.list <- list(Nomal = Normal.celltype, TypeIII = TypeIII.celltype)
#cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

# Define the cell labels to lift up, if the idents are overlaped, but can't cover each other
object.list <- list(Nomal = Normal.celltype, TypeIII = TypeIII.celltype)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

Normal.celltype <- liftCellChat(Normal.celltype , group.new = levels( cellchat@idents$joint) )

TypeIII.celltype <- liftCellChat(TypeIII.celltype , group.new = levels( cellchat@idents$joint) )
levels(Normal.celltype@idents)
levels(TypeIII.celltype@idents)



###############################################################
#Compare the total number of interactions and interaction strength
#To answer on question on whether the cell-cell communication is enhanced or not, 
#CellChat compares the the total number of interactions and interaction strength of 
#the inferred cell-cell communication networks from different biological conditions.
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

#To identify the interaction between which cell populations showing significant changes, 
#CellChat compares the number of interactions and interaction strength among different cell populations.

#The differential number of interactions or interaction strength in the cell-cell communication network 
#between two datasets can be visualized using circle plot, where red (or blue) colored edges represent 
#increased (or decreased) signaling in the second dataset compared to the first one.
par(mfrow = c(1,1), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

#We can also show differential number of interactions or interaction strength in a greater details using 
#a heatmap. The top colored bar plot represents the sum of column of values displayed in the heatmap 
#(incoming signaling). The right colored bar plot represents the sum of row of values (outgoing signaling). 
#In the colorbar, red (or blue) represents increased (or decreased) signaling in the second dataset compared 
#to the first one.
gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2
#To better control the node size and edge weights of the inferred networks across different datasets, we compute 
#the maximum number of cells per cell group and the maximum number of interactions (or interaction weights) 
#across all datasets.
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,1), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

#Comparing the outgoing and incoming interaction strength in 2D space allows ready identification of the cell 
#populations with significant changes in sending or receiving signals between different datasets.
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)

#Identify signaling groups based on their functional similarit
#library(umap)
#cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#cellchat <- netEmbedding(cellchat, type = "functional")
#cellchat <- netClustering(cellchat, type = "functional")
#netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)


#> Compute the distance of signaling networks between datasets 1 2
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

#Part III: Identify the upgulated and down-regulated signaling ligand-receptor pairs

# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "TypeIII"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "TypeIII",ligand.logFC = 0.2, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "Normal",ligand.logFC = -0.1, receptor.logFC = -0.1)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

#We then visualize the upgulated and down-regulated signaling ligand-receptor pairs using bubble plot or chord diagram
pairLR.use.up = net.up[, "interaction_name", drop = F]
#c(1,3,4,6,8,11,12,13,14,17,20,21,22), targets.use = c(2,5,7,9,10,15,16,18,19)
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c(1,4,7,8,11,12,15,17,18,19), targets.use = c(2,3,5,6,9,10,13,14,16), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c(2,3,5,6,9,10,13,14,16), targets.use = c(1,4,7,8,11,12,15,17,18,19), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2

#Visualize the upgulated and down-regulated signaling ligand-receptor pairs using Chord diagram
# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = c(1,4,7,8,11,12,15,17,18,19), targets.use = c(2,3,5,6,9,10,13,14,16),  slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = c(2,3,5,6,9,10,13,14,16), targets.use = c(1,4,7,8,11,12,15,17,18,19), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))


pathways.show <- c("MHC-I") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
# Chord diagram
pathways.show <- c("MHC-I") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}

pdf(file ="MHC-I.pdf", width = 10, height =8)
netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
dev.off()



#Part V: Compare the signaling gene expression distribution between different datasets
#We can plot the gene expression distribution of signaling genes related to L-R pairs or signaling pathway using a Seurat wrapper function plotGeneExpression.
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("wt", "ko")) # set factor level
plotGeneExpression(cellchat, signaling = "MHC-I", split.by = "datasets", colors.ggplot = T)
