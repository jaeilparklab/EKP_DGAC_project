Sys.setenv(LANG = "en_US")
library(Seurat)
devtools::install_github('satijalab/seurat-data')
library(SeuratData)
library(ggplot2)
library(sctransform)
library(dplyr)
library(patchwork)
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")

library(SeuratDisk)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DT)
library(biomaRt)
devtools::install_github('sunduanchen/Scissor')
library(Scissor)
memory.limit(512000)
options(future.globals.maxSize = 256000*1024^2)

# 1. prepare mouse gene symbol annotated TCGA-STAD RNA-seq gene matrix
# download TCGA-STAD gene matrix including 'Primary Tumor' and 'Solid Tissue Normal'
query <- GDCquery(project = "TCGA-STAD", data.category = "Transcriptome Profiling", data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts", legacy = FALSE, access= 'open', data.format='TXT', experimental.strategy='RNA-Seq', sample.type=c('Primary Tumor', 'Solid Tissue Normal'))
GDCdownload(query)
RNA_seq <- GDCprepare(query=query, save=TRUE, save.filename='/Users/gzou/OneDrive - Inside MD Anderson/Gengyi_MSO/Figures/Sissor/EKP_WT/TCGA_STAD_RNA_seq.rda')
GeneMatrix <- assay(RNA_seq)
write.csv(GeneMatrix, file="/Users/gzou/OneDrive - Inside MD Anderson/Gengyi_MSO/Figures/Sissor/EKP_WT/TCGA_STAD_GeneMatrix.csv", sep = ",", row.names=TRUE, col.names = TRUE)

# convert human ensembl id to mouse gene symbol
human <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", host='https://dec2021.archive.ensembl.org/')
mouse <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", host='https://dec2021.archive.ensembl.org/')
human_gene <- rownames(GeneMatrix)
mouse_gene <- getLDS(attributes = c("ensembl_gene_id_version"), filters = "ensembl_gene_id_version", values = human_gene, mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=TRUE)
head(mouse_gene)

mouse_gene <- mouse_gene[!(mouse_gene$MGI.symbol==""),]
mouse_gene <- mouse_gene[!duplicated(mouse_gene$Gene.stable.ID.version),]
mouse_gene <- mouse_gene[!duplicated(mouse_gene$MGI.symbol),]

# make the mouse gene symbol annotated gene matrix
class(GeneMatrix)
GeneMatrix_df<- as.data.frame(GeneMatrix)
class(GeneMatrix_df)

GeneMatrix_MGI <- GeneMatrix_df
GeneMatrix_MGI$genes <- mouse_gene$MGI.symbol[match(rownames(GeneMatrix_df), mouse_gene$Gene.stable.ID.version)]
GeneMatrix_MGI <- na.omit(GeneMatrix_MGI)
rownames(GeneMatrix_MGI) <- GeneMatrix_MGI$genes
GeneMatrix_MGI = GeneMatrix_MGI[,!(colnames(GeneMatrix_MGI) %in% 'genes')]
write.csv(GeneMatrix_MGI, file="/Users/gzou/OneDrive - Inside MD Anderson/Gengyi_MSO/Figures/Sissor/EKP_WT/GeneMatrix_MGI.csv", sep = ",", row.names=TRUE, col.names = TRUE)

# 2. prepare phenotype list
# download TCGA-STAD clinical index data
clinical_index <- GDCquery_clinic(project = "TCGA-STAD", type = "Clinical", save.csv = TRUE)

# Solid Tissue Normal
GeneMatrix_Normal <- subset(GeneMatrix, select = substr(colnames(GeneMatrix),14,15) %in% '11')
write.csv(GeneMatrix_Normal, file="/Users/gzou/OneDrive - Inside MD Anderson/Gengyi_MSO/Figures/Sissor/EKP_WT/GeneMatrix_Normal.csv", sep = ",", row.names=TRUE, col.names = TRUE)
Normal_phenotype <- as.data.frame(colnames(GeneMatrix_Normal))
Normal_phenotype$Disease <- 'Normal'
Normal_phenotype$Status <- 0
colnames(Normal_phenotype)<-  c("Patient", 'Disease', "Status")

# Carcinoma, diffuse type
diffuse <- clinical_index %>% filter(primary_diagnosis==c('Carcinoma, diffuse type'))
write.table(diffuse, file="/Users/gzou/OneDrive - Inside MD Anderson/Gengyi_MSO/Figures/Sissor/EKP_WT/diffuse_patients.csv", sep = ",", row.names=FALSE)
GeneMatrix_diffuse <- subset(GeneMatrix, select = substr(colnames(GeneMatrix),1,12) %in% diffuse$submitter_id)
GeneMatrix_diffuse <- subset(GeneMatrix_diffuse, select = substr(colnames(GeneMatrix_diffuse),14,15) %in% '01')
write.csv(GeneMatrix_diffuse, file="/Users/gzou/OneDrive - Inside MD Anderson/Gengyi_MSO/Figures/Sissor/EKP_WT/GeneMatrix_diffuse.csv", sep = ",", row.names=TRUE, col.names = TRUE)
diffuse_phenotype <- as.data.frame(colnames(GeneMatrix_diffuse))
diffuse_phenotype$Disease <- 'diffuse'
diffuse_phenotype$Status <- 1
colnames(diffuse_phenotype)<-  c("Patient", 'Disease', "Status")
Phenotype_normal_diffuse <- rbind(Normal_phenotype, diffuse_phenotype)
write.table(Phenotype_normal_diffuse, file="/Users/gzou/OneDrive - Inside MD Anderson/Gengyi_MSO/Figures/Sissor/EKP_WT/Phenotype_normal_diffuse.csv", sep = ",", row.names=FALSE)

# 3. match gene matrix and phenotype list
# only include patients that also exist in the phenotype list
GeneMatrix_normal_diffuse <- subset(GeneMatrix_MGI, select = colnames(GeneMatrix_MGI) %in% Phenotype_normal_diffuse$Patient)

# make gene matrix and phenotype list the same order
GeneMatrix_normal_diffuse <- GeneMatrix_normal_diffuse[ , order(match(colnames(GeneMatrix_normal_diffuse), Phenotype_normal_diffuse$Patient))]
all(colnames(GeneMatrix_normal_diffuse) == Phenotype_normal_diffuse$Patient)

# make the sample id annotated phenotype list
rownames(Phenotype_normal_diffuse) <- Phenotype_normal_diffuse$Patient
Phenotype_normal_diffuse_input = Phenotype_normal_diffuse[,!(colnames(Phenotype_normal_diffuse) %in% c('Patient','Disease'))]
Phenotype_normal_diffuse_input
names(Phenotype_normal_diffuse_input) = rownames(Phenotype_normal_diffuse)
Phenotype_normal_diffuse_input
table(Phenotype_normal_diffuse_input)
write.csv(Phenotype_normal_diffuse_input, file="/Users/gzou/OneDrive - Inside MD Anderson/Gengyi_MSO/Figures/Sissor/EKP_WT/Phenotype_normal_diffuse_input.csv", sep = ",", row.names=TRUE, col.names = TRUE)

# 4. prepare scRNA-seq data
# merge EKP and WT
EKP_data <- Read10X(data.dir = "/Users/gzou/OneDrive - Inside MD Anderson/Gengyi_MSO/MSO_matrix/EKP/filtered_feature_bc_matrix/")
EKP_data=EKP_data$`Gene Expression`
EKP <- CreateSeuratObject(counts = EKP_data, project = "EKP", min.cells = 3, min.features = 100)

WT_data <- Read10X(data.dir = "/Users/gzou/OneDrive - Inside MD Anderson/Gengyi_MSO/MSO_matrix/WT/filtered_feature_bc_matrix/")
WT_data=WT_data$`Gene Expression`
WT <- CreateSeuratObject(counts = WT_data, project = "WT", min.cells = 3, min.features = 100)

EKP$Batch<-"EKP"
WT$Batch<-"WT"
Merge_EKP_WT<-merge(x=EKP, y=c(WT))

Merge_EKP_WT[["percent.mt"]] <- PercentageFeatureSet(Merge_EKP_WT, pattern = "^mt-")
Merge_EKP_WT[["percent.rpl"]] <- PercentageFeatureSet(Merge_EKP_WT, pattern = "^Rpl")
Merge_EKP_WT[["percent.rps"]] <- PercentageFeatureSet(Merge_EKP_WT, pattern = "^Rps")

Merge_EKP_WT.P1 <- VlnPlot(Merge_EKP_WT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rpl","percent.rps"), ncol = 5)
Merge_EKP_WT.P1
Merge_EKP_WT.P2 <- FeatureScatter(Merge_EKP_WT, feature1 = "nCount_RNA", feature2 = "percent.mt")
Merge_EKP_WT.P3 <- FeatureScatter(Merge_EKP_WT, feature1 = "nCount_RNA", feature2 = "percent.rpl")
Merge_EKP_WT.P4 <- FeatureScatter(Merge_EKP_WT, feature1 = "nCount_RNA", feature2 = "percent.rps")
Merge_EKP_WT.P2+Merge_EKP_WT.P3+Merge_EKP_WT.P4
Merge_EKP_WT.P5 <- FeatureScatter(Merge_EKP_WT, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
Merge_EKP_WT.P5

Merge_EKP_WT <- NormalizeData(Merge_EKP_WT, normalization.method = "LogNormalize", scale.factor = 10000)
Merge_EKP_WT <- FindVariableFeatures(Merge_EKP_WT, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(Merge_EKP_WT)
Merge_EKP_WT <- ScaleData(Merge_EKP_WT, features = all.genes)
Merge_EKP_WT <- RunPCA(Merge_EKP_WT, verbose = FALSE)
Merge_EKP_WT <- FindNeighbors(Merge_EKP_WT, dims = 1:50, verbose = FALSE)
Merge_EKP_WT <- FindClusters(Merge_EKP_WT, resolution=0.5, verbose = FALSE)
Merge_EKP_WT <- RunUMAP(Merge_EKP_WT, dims = 1:50,verbose = FALSE)

Merge_EKP_WT_UMAP_p1<-DimPlot(Merge_EKP_WT, reduction = "umap",group.by = "Batch")
Merge_EKP_WT_UMAP_p1
Merge_EKP_WT_UMAP_p2<-DimPlot(Merge_EKP_WT, reduction = "umap",split.by = "Batch",label = TRUE,repel = TRUE,label.size = 13)+ NoLegend()
Merge_EKP_WT_UMAP_p2
Merge_EKP_WT_UMAP_p3<-DimPlot(Merge_EKP_WT, reduction = "umap", label = TRUE, repel = TRUE,label.size = 13)+ NoLegend()
Merge_EKP_WT_UMAP_p3

# 5. run Scissor
# EKP vs WT
phenotype <- Phenotype_normal_diffuse_input
tag <- c('Normal', 'diffuse')
GeneMatrix_normal_diffuse_input  <- as.matrix(GeneMatrix_normal_diffuse)
infos <- Scissor(bulk_dataset=GeneMatrix_normal_diffuse_input, sc_dataset=Merge_EKP_WT, phenotype=phenotype, tag = tag, alpha = 0.32, family = "binomial", Save_file = "normal_diffuse.RData")
Scissor_select <- rep(0, ncol(Merge_EKP_WT))
names(Scissor_select) <- colnames(Merge_EKP_WT)
Scissor_select[infos$Scissor_pos] <- 1
Scissor_select[infos$Scissor_neg] <- 2
sc_dataset <- AddMetaData(Merge_EKP_WT, metadata = Scissor_select, col.name = "scissor")
Scissor_EKP_WT_diffuse <- DimPlot(sc_dataset, reduction = 'umap', group.by = 'scissor', cols = c('grey','red','blue'), order = c(2,1))
Scissor_EKP_WT_diffuse
write.csv(sc_dataset[[]], file="/Users/gzou/OneDrive - Inside MD Anderson/Gengyi_MSO/Figures/Sissor/EKP_WT/scissor.csv", sep = ",", row.names=FALSE, col.names = FALSE)    # export meta.data of sc_dataset, you can check Scissor information
