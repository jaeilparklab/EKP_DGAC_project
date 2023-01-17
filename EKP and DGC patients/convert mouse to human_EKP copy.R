# convert mouse gene symbol of adata.var to human gene symbol, run below codes in R
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")
install.packages("RSQLite")
install.packages("cli")
library(biomaRt)
human <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mouse <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
GeneMatrix=read.csv(file="/Users/gzou/OneDrive - Inside MD Anderson/Gengyi_MSO/convert_to_human_gene/EKP_var.csv", header = TRUE, sep = ",", quote = "\"")
mouse_gene <- GeneMatrix[,1]
human_gene <- getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = mouse_gene, mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=TRUE)
head(human_gene)

human_gene <- human_gene[!(human_gene$HGNC.symbol==""),]    # remove rows with empty human gene
human_gene <- human_gene[!duplicated(human_gene$MGI.symbol),]    # remove rows with duplicated mouse gene
human_gene <- human_gene[!duplicated(human_gene$HGNC.symbol),]    # remove rows with duplicated human gene

class(GeneMatrix)
GeneMatrix_df<- as.data.frame(GeneMatrix)
class(GeneMatrix_df)

GeneMatrix_new <- GeneMatrix_df
GeneMatrix_new$genes <- human_gene$HGNC.symbol[match(GeneMatrix_df[,1], human_gene$MGI.symbol)]    # copy gene symbols to gene expression matrix, matching gene symbol
write.csv(GeneMatrix_new, file="/Users/gzou/OneDrive - Inside MD Anderson/Gengyi_MSO/convert_to_human_gene/EKP_GeneMatrix_new.csv", sep = ",", row.names=TRUE, col.names = TRUE)    # gene matrix with gene_symbol



