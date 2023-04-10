
  
#Single-Cell mRNAseq Analysis - HNSC data - Puram et al., 2017
#Data download from NCBI GEO - Single cell transcriptomes and metadata (sample/ cell type):
  
#Import HNSc scRNAseq data - Use read.delim2 to import cell type annotations as "strings", does not show up as NA if use read_tsv

setwd("C:/Users/jinsu/OneDrive/Desktop/Jinsu/Desktop Updated as of Mar 15, 2023/Research/Research with Dr. Bose/Dataset/single cell reference")

install.packages("Seurat")
install.packages("SeuratObject")
library(Seurat)
library(gdata)

library(EnhancedVolcano)

all_data <- read.delim2("HNSCC_all_data.txt")
#Clean Gene names - drop '' around gene name

colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HN26_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HN25_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC26_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HN28_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HN23_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC25_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC17_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC20_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC25_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC25_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC12_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC10_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC6_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC5_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC8_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC7_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC13_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC16_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC18_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC22_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC_17_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC_24_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC_28_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC24_", "", x)))




colnames(all_data)





colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P5", "LNM_YES_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P6", "LNM_YES_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P7", "LNM_YES_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P8", "LNM_NO_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P9", "LNM_NO_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P10", "LNM_YES_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P12", "LNM_NO_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P13", "LNM_YES_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P16", "LNM_NO_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P17", "LNM_NO_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P18", "LNM_YES_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P20", "LNM_YES_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P22", "LNM_NO_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P23", "LNM_YES_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P24", "LNM_YES_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P25", "LNM_YES_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P26", "LNM_YES_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P28", "LNM_YES_", x)))

colnames(all_data)
options(max.print = 5000)

all_data$X  <- as.list(sapply(all_data$X , function(x) gsub("\'", "", x)))
rownames(all_data) <- all_data$X
all_data$X <- NULL


index_TRUE <- startsWith(colnames(all_data), "LNM_YES_")
index_FALSE <- startsWith(colnames(all_data), "LNM_NO_")
LNM_PRESENT <-subset(all_data, select = ((startsWith(colnames(all_data), "LNM_YES_")== TRUE )))
LNM_ABSENT <-subset(all_data, select = ((startsWith(colnames(all_data), "LNM_NO_")== TRUE)))

#Create rna_meta data using the first 5 rows of the dataframe
rna_meta_PRESENT <- LNM_PRESENT[1:5,]
rna_meta_ABSENT <- LNM_ABSENT[1:5,]

rna_meta_PRESENT <- as.data.frame(t(rna_meta_PRESENT))
rna_meta_ABSENT <- as.data.frame(t(rna_meta_ABSENT))
rownames(rna_meta_PRESENT) <- gsub("_", "-", rownames(rna_meta_PRESENT))
rownames(rna_meta_ABSENT) <- gsub("_", "-", rownames(rna_meta_ABSENT))

index1 <- rna_meta_PRESENT$`non-cancer cell type` == 0
index2 <- rna_meta_ABSENT$`non-cancer cell type` == 0

rna_meta_PRESENT$`non-cancer cell type`[index1] <- "Cancer cell"
rna_meta_ABSENT$`non-cancer cell type`[index2] <- "Cancer cell"

index3 <- rna_meta_PRESENT$`non-cancer cell type` == "-Fibroblast"
index4 <- rna_meta_ABSENT$`non-cancer cell type` == "-Fibroblast"

rna_meta_PRESENT$`non-cancer cell type`[index3] <- "Fibroblast"
rna_meta_ABSENT$`non-cancer cell type`[index4] <- "Fibroblast"


#Create the counts matrix using the remaining data
rna_clean_PRESENT <- LNM_PRESENT[6:23691,]
rna_clean_ABSENT <- LNM_ABSENT[6:23691,]

names(rna_clean_PRESENT) <- gsub("_", "-", names(rna_clean_PRESENT))
names(rna_clean_ABSENT) <- gsub("_", "-", names(rna_clean_ABSENT))

#create seurat objects used for downstream analysis 
rna_seurat_LNM_POSITIVE <- CreateSeuratObject(rna_clean_PRESENT, project = "puram_data_LNM_POSITIVE")
rna_seurat_LNM_NEGATIVE <- CreateSeuratObject(rna_clean_ABSENT, project = "puram_data_LNM_ABSENT")

rna_seurat_LNM_POSITIVE <- AddMetaData(rna_seurat_LNM_POSITIVE, metadata = rna_meta_PRESENT)
rna_seurat_LNM_NEGATIVE <- AddMetaData(rna_seurat_LNM_NEGATIVE, metadata = rna_meta_ABSENT)



#Perfoming QC and selecting cells for further analysis

rna_seurat_LNM_POSITIVE[["percent.mt"]] <- PercentageFeatureSet(rna_seurat_LNM_POSITIVE, patter =
                                                                  "^MT-")
rna_seurat_LNM_NEGATIVE[["percent.mt"]] <- PercentageFeatureSet(rna_seurat_LNM_NEGATIVE, patter =
                                                                  "^MT-")

VlnPlot(rna_seurat_LNM_POSITIVE, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(rna_seurat_LNM_NEGATIVE, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


plot1 <- FeatureScatter(rna_seurat_LNM_POSITIVE, feature1 = "nCount_RNA", feature2= "percent.mt")
plot2 <- FeatureScatter(rna_seurat_LNM_POSITIVE, feature1 = "nCount_RNA", feature2= "nFeature_RNA")
plot3 <- FeatureScatter(rna_seurat_LNM_NEGATIVE, feature1 = "nCount_RNA", feature2= "percent.mt")
plot4 <- FeatureScatter(rna_seurat_LNM_NEGATIVE, feature1 = "nCount_RNA", feature2= "nFeature_RNA")


rna_seurat_LNM_POSITIVE <- NormalizeData(rna_seurat_LNM_POSITIVE, 
                                         normalization.method = "LogNormalize", 
                                         scale.factor = 10000)
rna_seurat_LNM_NEGATIVE <- NormalizeData(rna_seurat_LNM_NEGATIVE, 
                                        normalization.method = "LogNormalize", 
                                        scale.factor = 10000)

rna_seurat_LNM_POSITIVE <- FindVariableFeatures(rna_seurat_LNM_POSITIVE, selection.method = "vst",
                                                nfeatures = 9000)
 
rna_seurat_LNM_NEGATIVE <- FindVariableFeatures(rna_seurat_LNM_NEGATIVE, selection.method = "vst",
                                                nfeatures = 9000)   

top10_P <- head(VariableFeatures(rna_seurat_LNM_POSITIVE))
top10_A <- head(VariableFeatures(rna_seurat_LNM_NEGATIVE))

plot5 <- VariableFeaturePlot(rna_seurat_LNM_POSITIVE)
plot6 <- VariableFeaturePlot(rna_seurat_LNM_NEGATIVE)


all.genes_P <- rownames(rna_seurat_LNM_POSITIVE)
all.genes_N <- rownames(rna_seurat_LNM_NEGATIVE)

rna_seurat_LNM_POSITIVE <- ScaleData(rna_seurat_LNM_POSITIVE, features = all.genes_P)
rna_seurat_LNM_NEGATIVE <- ScaleData(rna_seurat_LNM_NEGATIVE, features = all.genes_N)


rna_seurat_LNM_POSITIVE <- RunPCA(rna_seurat_LNM_POSITIVE, 
                                  features = VariableFeatures(object =rna_seurat_LNM_POSITIVE))
rna_seurat_LNM_NEGATIVE <- RunPCA(rna_seurat_LNM_NEGATIVE, 
                                  features = VariableFeatures(object =rna_seurat_LNM_NEGATIVE))


print(rna_seurat_LNM_POSITIVE[["pca"]], dims = 1:5, nfeatures = 5)
print(rna_seurat_LNM_NEGATIVE[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(rna_seurat_LNM_POSITIVE, dims = 1:2, reduction = "pca")
DimPlot(rna_seurat_LNM_POSITIVE, reduction = "pca")
DimHeatmap(rna_seurat_LNM_POSITIVE, dims = 1:15, cells = 500, balanced =  TRUE)

rna_seurat_LNM_POSITIVE <- JackStraw(rna_seurat_LNM_POSITIVE, num.replicate = 100)
rna_seurat_LNM_POSITIVE <- ScoreJackStraw(rna_seurat_LNM_POSITIVE, dims = 1:20)
JackStrawPlot(rna_seurat_LNM_POSITIVE, dims= 1:15)
ElbowPlot(rna_seurat_LNM_POSITIVE)

rna_seurat_LNM_POSITIVE <- FindNeighbors(rna_seurat_LNM_POSITIVE, dims = 1:15)
rna_seurat_LNM_POSITIVE <- FindClusters(rna_seurat_LNM_POSITIVE, resolution = 0.5)

head(Idents(rna_seurat_LNM_POSITIVE), 5)
rna_seurat_LNM_POSITIVE <- RunUMAP(rna_seurat_LNM_POSITIVE, dims = 1:15)
rna_seurat_LNM_POSITIVE <- RunTSNE(rna_seurat_LNM_POSITIVE, dims = 1:15)
DimPlot(rna_seurat_LNM_POSITIVE, reduction = "umap", group.by = 'non.cancer.cell.type')
DimPlot(rna_seurat_LNM_POSITIVE, reduction = "tsne", group.by = 'non.cancer.cell.type')


rna_seurat_LNM_NEGATIVE <- FindNeighbors(rna_seurat_LNM_NEGATIVE, dims = 1:15)
rna_seurat_LNM_NEGATIVE <- FindClusters(rna_seurat_LNM_NEGATIVE, resolution = 0.3)

head(Idents(rna_seurat_LNM_NEGATIVE), 5)
rna_seurat_LNM_NEGATIVE <- RunUMAP(rna_seurat_LNM_NEGATIVE, dims = 1:15)
rna_seurat_LNM_NEGATIVE <- RunTSNE(rna_seurat_LNM_NEGATIVE, dims = 1:15)
DimPlot(rna_seurat_LNM_NEGATIVE, reduction = "umap", group.by = 'non.cancer.cell.type')
DimPlot(rna_seurat_LNM_NEGATIVE, reduction = "tsne", group.by = 'non.cancer.cell.type')

#DEGs for Cancer cell type
Idents(rna_seurat_LNM_POSITIVE) <- 'non.cancer.cell.type'
Idents(rna_seurat_LNM_NEGATIVE) <- 'non.cancer.cell.type'

Cancer_Subsetted_LNM_PRESENT <- subset(x = rna_seurat_LNM_POSITIVE, idents = "Cancer cell")
Cancer_Subsetted_LNM_ABSENT <- subset(x = rna_seurat_LNM_NEGATIVE, idents = "Cancer cell")


Idents(object = Cancer_Subsetted_LNM_PRESENT) <- "Metatasis"
Idents(object = Cancer_Subsetted_LNM_ABSENT) <- "Non-Metatasis"


Merged_Cancer_celltype <- merge(x = Cancer_Subsetted_LNM_PRESENT, y = Cancer_Subsetted_LNM_ABSENT, project = "merged_SCC_seurat")
DEGs_inCancer_LNMPOS <- FindMarkers(Merged_Cancer_celltype, ident.1 = "Metatasis", ident.2 = "Non-Metatasis", recorrect_umi = FALSE)

#DEGs for T cell type
T_Subsetted_LNM_PRESENT <- subset(x = rna_seurat_LNM_POSITIVE, idents = "T cell")
T_Subsetted_LNM_ABSENT <- subset(x = rna_seurat_LNM_NEGATIVE, idents = "T cell")


Idents(object = T_Subsetted_LNM_PRESENT) <- "Metatasis"
Idents(object = T_Subsetted_LNM_ABSENT) <- "Non-Metatasis"


Merged_T_cell_celltype <- merge(x = T_Subsetted_LNM_PRESENT, y = T_Subsetted_LNM_ABSENT, project = "merged_SCC_seurat")
DEGs_in_Tcells_LNMPOS <- FindMarkers(Merged_T_cell_celltype, ident.1 = "Metatasis", ident.2 = "Non-Metatasis", recorrect_umi = FALSE)

#DEGs for Fibroblast cell type
Fibroblast_Subsetted_LNM_PRESENT <- subset(x = rna_seurat_LNM_POSITIVE, idents = "Fibroblast")
Fibroblast_Subsetted_LNM_ABSENT <- subset(x = rna_seurat_LNM_NEGATIVE, idents = "Fibroblast")

Idents(object = Fibroblast_Subsetted_LNM_PRESENT) <- "Metatasis"
Idents(object = Fibroblast_Subsetted_LNM_ABSENT) <- "Non-Metatasis"

Merged_Fibroblast_celltype <- merge(x = Fibroblast_Subsetted_LNM_PRESENT, y = Fibroblast_Subsetted_LNM_ABSENT, project = "merged_SCC_seurat")
DEGs_in_Fibroblast_LNMPOS <- FindMarkers(Merged_Fibroblast_celltype, ident.1 = "Metatasis", ident.2 = "Non-Metatasis", recorrect_umi = FALSE)

Bcell_Subsetted_LNM_PRESENT <- subset(x = rna_seurat_LNM_POSITIVE, idents = "B cell")
Bcell_Subsetted_LNM_ABSENT <- subset(x = rna_seurat_LNM_NEGATIVE, idents = "B cell")

Idents(object = Bcell_Subsetted_LNM_PRESENT) <- "Metatasis"
Idents(object = Bcell_Subsetted_LNM_ABSENT) <- "Non-Metatasis"

Merged_Bcell_celltype <- merge(x = Bcell_Subsetted_LNM_PRESENT, y = Bcell_Subsetted_LNM_ABSENT, project = "merged_SCC_seurat")
DEGs_in_Bcell_LNMPOS <- FindMarkers(Merged_Bcell_celltype, ident.1 = "Metatasis", ident.2 = "Non-Metatasis", recorrect_umi = FALSE)

Macrophage_Subsetted_LNM_PRESENT <- subset(x = rna_seurat_LNM_POSITIVE, idents = "Macrophage")
Macrophage_Subsetted_LNM_ABSENT <- subset(x = rna_seurat_LNM_NEGATIVE, idents = "Macrophage")

Idents(object = Macrophage_Subsetted_LNM_PRESENT) <- "Metatasis"
Idents(object = Macrophage_Subsetted_LNM_ABSENT) <- "Non-Metatasis"

Merged_Macrophage_celltype <- merge(x = Macrophage_Subsetted_LNM_PRESENT, y = Macrophage_Subsetted_LNM_ABSENT, project = "merged_SCC_seurat")
DEGs_in_Macrophage_LNMPOS <- FindMarkers(Merged_Macrophage_celltype, ident.1 = "Metatasis", ident.2 = "Non-Metatasis", recorrect_umi = FALSE)

Endothelial_Subsetted_LNM_PRESENT <- subset(x = rna_seurat_LNM_POSITIVE, idents = "Endothelial")
Endothelial_Subsetted_LNM_ABSENT <- subset(x = rna_seurat_LNM_NEGATIVE, idents = "Endothelial")

Idents(object = Endothelial_Subsetted_LNM_PRESENT) <- "Metatasis"
Idents(object = Endothelial_Subsetted_LNM_ABSENT) <- "Non-Metatasis"

Merged_Endothelial_celltype <- merge(x = Endothelial_Subsetted_LNM_PRESENT, y = Endothelial_Subsetted_LNM_ABSENT, project = "merged_SCC_seurat")
DEGs_in_Endothelial_LNMPOS <- FindMarkers(Merged_Endothelial_celltype, ident.1 = "Metatasis", ident.2 = "Non-Metatasis", recorrect_umi = FALSE)


Mast_Subsetted_LNM_PRESENT <- subset(x = rna_seurat_LNM_POSITIVE, idents = "Mast")
Mast_Subsetted_LNM_ABSENT <- subset(x = rna_seurat_LNM_NEGATIVE, idents = "Mast")

Idents(object = Mast_Subsetted_LNM_PRESENT) <- "Metatasis"
Idents(object = Mast_Subsetted_LNM_ABSENT) <- "Non-Metatasis"

Merged_Mast_celltype <- merge(x = Mast_Subsetted_LNM_PRESENT, y = Mast_Subsetted_LNM_ABSENT, project = "merged_SCC_seurat")
DEGs_in_Mast_LNMPOS <- FindMarkers(Merged_Mast_celltype, ident.1 = "Metatasis", ident.2 = "Non-Metatasis", recorrect_umi = FALSE)


View(DEGs_inCancer_LNMPOS)
View(DEGs_in_Tcells_LNMPOS)
View(DEGs_in_Fibroblast_LNMPOS)


p_val_Cancer <- subset(DEGs_inCancer_LNMPOS, p_val_adj < 0.01)
logFold_greaterthan0.58_cancer <- subset(p_val_Cancer,  avg_log2FC > 0.25 )
logFold_lessthan0.58_cancer <- subset(p_val_Cancer,  avg_log2FC < -0.25)
Filtered_Cancer <- combine(logFold_greaterthan0.58_cancer, logFold_lessthan0.58_cancer)

p_val_T <- subset(DEGs_in_Tcells_LNMPOS, p_val_adj < 0.05)
logFold_greaterthan0.58_T <- subset(p_val_T,  avg_log2FC > 0.25 )
logFold_lessthan0.58_T <- subset(p_val_T,  avg_log2FC < -0.25)
Filtered_T <- combine(logFold_greaterthan0.58_T, logFold_lessthan0.58_T)

p_val_Fibroblast <- subset(DEGs_in_Fibroblast_LNMPOS, p_val_adj < 0.05)
logFold_greaterthan0.58_Fibroblast <- subset(p_val_Fibroblast,  avg_log2FC > 0.25 )
logFold_lessthan0.58_Fibroblast <- subset(p_val_Fibroblast,  avg_log2FC < -0.25)
Filtered_Fibroblast <- combine(logFold_greaterthan0.58_Fibroblast, logFold_lessthan0.58_Fibroblast)

p_val_B <- subset(DEGs_in_Bcell_LNMPOS, p_val_adj < 0.05)
logFold_greaterthan0.58_B <- subset(p_val_B,  avg_log2FC > 0.25 )
logFold_lessthan0.58_B <- subset(p_val_B,  avg_log2FC < -0.25)
Filtered_B <- combine(logFold_greaterthan0.58_B, logFold_lessthan0.58_B)

p_val_Macrophage <- subset(DEGs_in_Macrophage_LNMPOS, p_val_adj < 0.05)
logFold_greaterthan0.58_Macrophage <- subset(p_val_Macrophage,  avg_log2FC > 0.25 )
logFold_lessthan0.58_Macrophage <- subset(p_val_Macrophage,  avg_log2FC < -0.25)
Filtered_Macrophage <- combine(logFold_greaterthan0.58_Macrophage, logFold_lessthan0.58_Macrophage)

p_val_Mast <- subset(DEGs_in_Mast_LNMPOS, p_val_adj < 0.05)
logFold_greaterthan0.58_Mast <- subset(p_val_Mast,  avg_log2FC > 0.25 )
logFold_lessthan0.58_Mast <- subset(p_val_Mast,  avg_log2FC < -0.25)
Filtered_Mast <- combine(logFold_greaterthan0.58_Mast, logFold_lessthan0.58_Mast)

write.csv(Filtered_Cancer, file = "C:/Users/jinsu/OneDrive/Desktop/Jinsu/CancerDEG.csv")

saveRDS(rna_seurat_LNM_POSITIVE, file = "rna_seurat_LNM_POSITIVE")
saveRDS(rna_seurat_LNM_NEGATIVE, file = "rna_seurat_LNM_NEGATIVE")

BiocManager::install("DESeq2")
library(DESeq2)

BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

BiocManager::install("goseq", force =TRUE)

BiocManager::install("org.Hs.eg.db", force =TRUE)
BiocManager::install("dplyr",force =TRUE )
BiocManager::install("fgsea", force =TRUE)
BiocManager::install("ggplot2", force =TRUE)
BiocManager::install("pheatmap", force =TRUE)
BiocManager::install("RColorBrewer", force =TRUE)
BiocManager::install("gridExtra", force =TRUE)
BiocManager::install("rmarkdown", force =TRUE)


library(goseq)
library(org.Hs.eg.db)
library(dplyr)
library(fgsea)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(gridExtra)
library(rmarkdown)


BiocManager::install("GenomicRanges", force = TRUE)

library(Seurat)

install.packages("GRanges")
library(GenomicRanges)
#1
markers <- FindMarkers(pbmc, ident.1 = 'NK', only.pos = T)


head(markers, n = 10)

pbmc <- readRDS("pbmc3k_final.rds")


supportedOrganisms()[supportedOrganisms()$Genome=="hg38",]
supportedOrganisms()[supportedOrganisms()$Genome=="hg19",]

all.genes <- rownames(rna_seurat_LNM_POSITIVE)

de_genes <- rownames(Filtered_Cancer)

gene_vector <- as.integer(all.genes%in%de_genes)

names(gene_vector) <- all.genes

head(gene_vector)

pwf <- nullp(gene_vector, "hg19", "geneSymbol")

go_wall <- goseq(pwf, "hg19", "geneSymbol")

head(go_wall, 20)
paged_table(head(go_wall, 20))

go_wall %>% 
  top_n(10, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=term, 
             colour=over_represented_pvalue, 
             size=numDEInCat), labSize = 4.0) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count")



ranks  <- Filtered_Cancer$avg_log2FC
names(ranks) <- rownames(Filtered_Cancer)

barplot(sort(ranks, decreasing = T))

pathways_hallmark <- gmtPathways("h.all.v7.2.symbols.gmt")


pathways_hallmark %>% 
  head() %>% 
  lapply(head)

fgsea_res <- fgsea(pathways=pathways_hallmark, stats=ranks, nperm=1000)

# check results
head(fgsea_res[order(padj, -abs(NES)), ], n=10)

# plot
fgsea_res_tidy <- fgsea_res %>%  as_tibble() %>% arrange(desc(NES))

p <- ggplot(fgsea_res_tidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

plot(p)  

plotEnrichment(pathways_hallmark[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]], ranks)
plotEnrichment(pathways_hallmark[["HALLMARK_ANGIOGENESIS"]], ranks)


DEGs_NO_FILTER <- EnhancedVolcano(DEGs_inCancer_LNMPOS , 
                                  lab = rownames(DEGs_inCancer_LNMPOS),
                                  x ="avg_log2FC", 
                                  y ="p_val_adj",
                                  title = 'DEGs in Cancer cell type of LNM+ before Filtering ',
                                  pCutoff = 0.01,
                                  FCcutoff = 0.25,
                                  pointSize = 1.5,
                                  labSize = 4.0)

DEGs_Filter <- EnhancedVolcano(Filtered_Cancer , 
                                  lab = rownames(Filtered_Cancer),
                                  x ="avg_log2FC", 
                                  y ="p_val_adj",
                                  title = 'DEGs in Cancer cell type of LNM+ after Filtering ',
                                  pCutoff = 0.01,
                                  FCcutoff = 0.25,
                                  pointSize = 1.5,
                                  labSize = 4.0)


DEGs_Filter <- EnhancedVolcano(DEG , 
                               lab = rownames(DEGs_in_Mast_LNMPOS),
                               x ="avg_log2FC", 
                               y ="p_val_adj",
                               title = 'DEGs in Mast cell type of LNM+ after Filtering ',
                               pCutoff = 0.01, 
                               FCcutoff = 0.25,
                               pointSize = 1.5,
                               labSize = 4.0)

DEGs_Filter


















#Perform basic data validation - examine basic characteristics of data
counts_per_cell_P <- Matrix::colSums(rna_seurat_LNM_POSITIVE)
counts_per_gene_P <- Matrix::rowSums(rna_seurat_LNM_POSITIVE)
genes_per_cell_P <- Matrix::colSums(rna_seurat_LNM_POSITIVE)
cells_per_gene_P <- Matrix::rowSums(rna_seurat_LNM_POSITIVE)

counts_per_cell_A <- Matrix::colSums(rna_seurat_LNM_NEGATIVE)
counts_per_gene_A <- Matrix::rowSums(rna_seurat_LNM_NEGATIVE)
genes_per_cell_A <- Matrix::colSums(rna_seurat_LNM_NEGATIVE)
cells_per_gene_A <- Matrix::rowSums(rna_seurat_LNM_NEGATIVE)


#Graph the counts per cell, counts per gene, cells per gene, etc 
hist(log10(counts_per_cell+1),main='counts per cell',col='wheat')
hist(log10(genes_per_cell+1), main='genes per cell', col='wheat')
plot(counts_per_cell, genes_per_cell, log='xy', col='wheat')
title('counts vs genes per cell')
hist(log10(counts_per_gene+1), main='counts per gene', col='wheat')
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')



#Normalize the RNA data in the seurat object:
rna_seurat_LNM_POSITIVE <- NormalizeData(rna_seurat_LNM_POSITIVE)
rna_seurat_LNM_NEGATIVE <- NormalizeData(rna_seurat_LNM_NEGATIVE)

rna_seurat_LNM_POSITIVE <- FindVariableFeatures(rna_seurat_LNM_POSITIVE)
rna_seurat_LNM_NEGATIVE <- FindVariableFeatures(rna_seurat_LNM_NEGATIVE)


top10_P <- head(VariableFeatures(rna_seurat_LNM_POSITIVE), 10)
top10_A <- head(VariableFeatures(rna_seurat_LNM_NEGATIVE), 10)

plot1_P <- VariableFeaturePlot(rna_seurat_LNM_POSITIVE)
LabelPoints(plot = plot1_P, points = top10_P, repel = TRUE, xnudge = 0, ynudge = 0)

plot1_A <- VariableFeaturePlot(rna_seurat_LNM_NEGATIVE)
LabelPoints(plot = plot1_A, points = top10_A, repel = TRUE, xnudge = 0, ynudge = 0)


rna_seurat_LNM_POSITIVE <- ScaleData(rna_seurat_LNM_POSITIVE)
rna_seurat_LNM_NEGATIVE <- ScaleData(rna_seurat_LNM_NEGATIVE)


rna_seurat_LNM_POSITIVE <- RunPCA(rna_seurat_LNM_POSITIVE, verbose = FALSE)
rna_seurat_LNM_NEGATIVE <- RunPCA(rna_seurat_LNM_NEGATIVE, verbose = FALSE)


DimHeatmap(rna_seurat_LNM_POSITIVE, dims = 1:6, nfeatures = 20, cells = 500, balanced = T)
ElbowPlot(rna_seurat_LNM_POSITIVE)

DimHeatmap(rna_seurat_LNM_NEGATIVE, dims = 1:6, nfeatures = 20, cells = 500, balanced = T)
ElbowPlot(rna_seurat_LNM_NEGATIVE)

rna_seurat_LNM_POSITIVE <- FindNeighbors(rna_seurat_LNM_POSITIVE, dims = 1:10)
rna_seurat_LNM_NEGATIVE <- FindNeighbors(rna_seurat_LNM_NEGATIVE, dims = 1:10)

rna_seurat_LNM_POSITIVE <- FindClusters(rna_seurat_LNM_POSITIVE, resolution = 0.8, verbose = FALSE)
rna_seurat_LNM_NEGATIVE <- FindClusters(rna_seurat_LNM_NEGATIVE, resolution = 0.8, verbose = FALSE)

rna_seurat_LNM_POSITIVE <- RunUMAP(rna_seurat_LNM_POSITIVE, dims = 1:10, verbose = FALSE)
rna_seurat_LNM_NEGATIVE <- RunUMAP(rna_seurat_LNM_NEGATIVE, dims = 1:10, verbose = FALSE)

DimPlot(rna_seurat_LNM_POSITIVE, reduction = "umap", group.by = "non.cancer.cell.type", label = TRUE, repel =TRUE, DarkTheme()) 
DimPlot(rna_seurat_LNM_NEGATIVE, reduction = "umap", group.by = "non.cancer.cell.type", label = TRUE, repel =TRUE) 

rna_seurat_LNM_POSITIVE <- RunTSNE(rna_seurat_LNM_POSITIVE, dims = 1:10)
rna_seurat_LNM_NEGATIVE <- RunTSNE(rna_seurat_LNM_NEGATIVE, dims = 1:10)

plot_for_LNM_POS <- DimPlot(rna_seurat_LNM_POSITIVE, reduction = "tsne", group.by = "non.cancer.cell.type", label = TRUE, repel =TRUE)
plot_for_LNM_NEG <- DimPlot(rna_seurat_LNM_NEGATIVE, reduction = "tsne", group.by = "non.cancer.cell.type", label = TRUE, repel =TRUE)

plot_for_LNM_POS + plot_for_LNM_NEG


#DEGs for Cancer cell type
Idents(rna_seurat_LNM_POSITIVE) <- 'non.cancer.cell.type'
Idents(rna_seurat_LNM_NEGATIVE) <- 'non.cancer.cell.type'

Cancer_Subsetted_LNM_PRESENT <- subset(x = rna_seurat_LNM_POSITIVE, idents = "Cancer cell")
Cancer_Subsetted_LNM_ABSENT <- subset(x = rna_seurat_LNM_NEGATIVE, idents = "Cancer cell")

top10_P <- head(VariableFeatures(Cancer_Subsetted_LNM_PRESENT), 10)
top10_A <- head(VariableFeatures(Cancer_Subsetted_LNM_ABSENT), 10)


Idents(object = Cancer_Subsetted_LNM_PRESENT) <- "Metatasis"
Idents(object = Cancer_Subsetted_LNM_ABSENT) <- "Non-Metatasis"


Merged_Cancer_celltype <- merge(x = Cancer_Subsetted_LNM_PRESENT, y = Cancer_Subsetted_LNM_ABSENT, project = "merged_SCC_seurat")
DEGs_inCancer_LNMPOS <- FindMarkers(Merged_Cancer_celltype, ident.1 = "Metatasis", ident.2 = "Non-Metatasis", recorrect_umi = FALSE)

#DEGs for T cell type
T_Subsetted_LNM_PRESENT <- subset(x = rna_seurat_LNM_POSITIVE, idents = "T cell")
T_Subsetted_LNM_ABSENT <- subset(x = rna_seurat_LNM_NEGATIVE, idents = "T cell")


Idents(object = T_Subsetted_LNM_PRESENT) <- "Metatasis"
Idents(object = T_Subsetted_LNM_ABSENT) <- "Non-Metatasis"


Merged_T_cell_celltype <- merge(x = T_Subsetted_LNM_PRESENT, y = T_Subsetted_LNM_ABSENT, project = "merged_SCC_seurat")
DEGs_in_Tcells_LNMPOS <- FindMarkers(Merged_T_cell_celltype, ident.1 = "Metatasis", ident.2 = "Non-Metatasis", recorrect_umi = FALSE)

#DEGs for Fibroblast cell type
Fibroblast_Subsetted_LNM_PRESENT <- subset(x = rna_seurat_LNM_POSITIVE, idents = "Fibroblast")
Fibroblast_Subsetted_LNM_ABSENT <- subset(x = rna_seurat_LNM_NEGATIVE, idents = "Fibroblast")

Idents(object = Fibroblast_Subsetted_LNM_PRESENT) <- "Metatasis"
Idents(object = Fibroblast_Subsetted_LNM_ABSENT) <- "Non-Metatasis"

Merged_Fibroblast_celltype <- merge(x = Fibroblast_Subsetted_LNM_PRESENT, y = Fibroblast_Subsetted_LNM_ABSENT, project = "merged_SCC_seurat")
DEGs_in_Fibroblast_LNMPOS <- FindMarkers(Merged_Fibroblast_celltype, ident.1 = "Metatasis", ident.2 = "Non-Metatasis", recorrect_umi = FALSE)

View(DEGs_inCancer_LNMPOS)
View(DEGs_in_Tcells_LNMPOS)
View(DEGs_in_Fibroblast_LNMPOS)

p_val_Cancer <- subset(DEGs_inCancer_LNMPOS, p_val_adj < 0.001)
logFold_greaterthan0.58_cancer <- subset(p_val_Cancer,  avg_log2FC > 0.58 )
logFold_lessthan0.58_cancer <- subset(p_val_Cancer,  avg_log2FC < -0.58)
Filtered_Cancer <- combine(logFold_greaterthan0.58_cancer, logFold_lessthan0.58_cancer)

p_val_T <- subset(DEGs_in_Tcells_LNMPOS, p_val_adj < 0.001)
logFold_greaterthan0.58_T <- subset(p_val_T,  avg_log2FC > 0.58 )
logFold_lessthan0.58_T <- subset(p_val_T,  avg_log2FC < -0.58)
Filtered_T <- combine(logFold_greaterthan0.58_T, logFold_lessthan0.58_T)

p_val_Fibroblast <- subset(DEGs_in_Fibroblast_LNMPOS, p_val_adj < 0.001)
logFold_greaterthan0.58_Fibroblast <- subset(p_val_Fibroblast,  avg_log2FC > 0.58 )
logFold_lessthan0.58_Fibroblast <- subset(p_val_Fibroblast,  avg_log2FC < -0.58)
Filtered_Fibroblast <- combine(logFold_greaterthan0.58_Fibroblast, logFold_lessthan0.58_Fibroblast)
















#####
table(rna_seurat@meta.data$seurat_clusters)
DimPlot(rna_seurat, label.size = 4, repel = TRUE, label = TRUE)

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

rna_seurat <- CellCycleScoring(rna_seurat, s.features = s.genes, g2m.features = g2m.genes)
table(rna_seurat[[]]$Phase)

#FeaturePlot(rna_seurat, features = c("LILRA4", "TPM2", "PPBP", "GP1BB"))

rna_seurat <- RunUMAP(rna_seurat, dims = 1:10)
DimPlot(rna_seurat, reduction = "umap", group.by = "non.cancer.cell.type", label = TRUE, repel =TRUE) 

rna_seurat <- RunTSNE(rna_seurat, dims = 1:10)
DimPlot(rna_seurat, reduction = "tsne", group.by = "non.cancer.cell.type", label = TRUE, repel =TRUE)

FeatureScatter(object = rna_seurat, feature1 = "RNF126", feature2 = "PC_1")

VlnPlot(object = rna_seurat, features = c("MIER2"))












all_data <- read.delim2("HNSCC_all_data.txt")

#Clean Gene names - drop '' around gene name
all_data$X  <- as.list(sapply(all_data$X , function(x) gsub("\'", "", x)))
rownames(all_data) <- all_data$X
all_data$X <- NULL
#Create rna_meta data using the first 5 rows of the dataframe
rna_meta <- all_data[1:5,]
rna_meta <- as.data.frame(t(rna_meta))
rownames(rna_meta) <- gsub("_", "-", rownames(rna_meta))


index <- rna_meta$`non-cancer cell type` == "-Fibroblast"
rna_meta$`non-cancer cell type`[index] <- "Fibroblast"

index <- rna_meta$`non-cancer cell type` == 0
rna_meta$`non-cancer cell type`[index] <- "Cancer cell"
#Create the counts matrix using the remaining data
rna_clean <- all_data[6:23691,]
names(rna_clean) <- gsub("_", "-", names(rna_clean))
#create seurat objects used for downstream analysis 
rna_seurat <- CreateSeuratObject(rna_clean, project = "puram_data")
rna_seurat <- AddMetaData(rna_seurat, metadata = rna_meta)

#Perform basic data validation - examine basic characteristics of data
counts_per_cell <- Matrix::colSums(rna_seurat)
counts_per_gene <- Matrix::rowSums(rna_seurat)
genes_per_cell <- Matrix::colSums(rna_seurat)
cells_per_gene <- Matrix::rowSums(rna_seurat)
#Graph the counts per cell, counts per gene, cells per gene, etc 
hist(log10(counts_per_cell+1),main='counts per cell',col='wheat')
hist(log10(genes_per_cell+1), main='genes per cell', col='wheat')
plot(counts_per_cell, genes_per_cell, log='xy', col='wheat')
title('counts vs genes per cell')
hist(log10(counts_per_gene+1), main='counts per gene', col='wheat')
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')

#Normalize the RNA data in the seurat object:
rna_seurat <- NormalizeData(rna_seurat)
rna_seurat <- FindVariableFeatures(rna_seurat)


top10 <- head(VariableFeatures(rna_seurat), 10)
top10

plot1 <- VariableFeaturePlot(rna_seurat)
LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)

rna_seurat <- ScaleData(rna_seurat)

rna_seurat <- RunPCA(rna_seurat, verbose = FALSE)
DimHeatmap(rna_seurat, dims = 1:6, nfeatures = 20, cells = 500, balanced = T)
ElbowPlot(rna_seurat)

rna_seurat <- FindNeighbors(rna_seurat, dims = 1:10)

rna_seurat <- FindClusters(rna_seurat, resolution = 0.8, verbose = FALSE)

rna_seurat <- RunUMAP(rna_seurat, dims = 1:10, verbose = FALSE)
table(rna_seurat@meta.data$seurat_clusters)
DimPlot(rna_seurat, label.size = 4, repel = TRUE, label = TRUE)

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

rna_seurat <- CellCycleScoring(rna_seurat, s.features = s.genes, g2m.features = g2m.genes)
table(rna_seurat[[]]$Phase)

#FeaturePlot(rna_seurat, features = c("LILRA4", "TPM2", "PPBP", "GP1BB"))

rna_seurat <- RunUMAP(rna_seurat, dims = 1:10)
DimPlot(rna_seurat, reduction = "umap")

rna_seurat <- RunTSNE(rna_seurat, dims = 1:10)
DimPlot(rna_seurat, reduction = "tsne")


Plot1 <- DimPlot(rna_seurat, group.by = 'Phase')
Plot2 <- DimPlot(rna_seurat, group.by = 'non.cancer.cell.type')


Plot1 + Plot2

FeatureScatter(object = rna_seurat, feature1 = "RNF126", feature2 = "PC_1")

VlnPlot(object = rna_seurat, features = c("MIER2"))




Idents(object = rna_seurat) <- "non.cancer.cell.type"

Fibrobalst <-subset(x = rna_seurat, idents = "Fibroblast")

DimPlot(Fibrobalst, group.by = 'Phase')




























#Clean Gene names - drop '' around gene name & only patient number available
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HN26_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HN25_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC26_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HN28_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HN23_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC25_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC17_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC20_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC25_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC25_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC12_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC10_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC6_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC5_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC8_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC7_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC13_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC16_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC18_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC22_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC_17_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC_24_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC_28_", "", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("HNSCC24_", "", x)))




colnames(all_data)





colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P5", "LNM_YES_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P6", "LNM_YES_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P7", "LNM_YES_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P8", "LNM_NO_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P9", "LNM_NO_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P10", "LNM_YES_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P12", "LNM_NO_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P13", "LNM_YES_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P16", "LNM_NO_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P17", "LNM_NO_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P18", "LNM_YES_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P20", "LNM_YES_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P22", "LNM_NO_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P23", "LNM_YES_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P24", "LNM_YES_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P25", "LNM_YES_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P26", "LNM_YES_", x)))
colnames(all_data) <- as.list(sapply(colnames(all_data), function(x) sub("P28", "LNM_YES_", x)))

colnames(all_data)
options(max.print = 5000)


index_TRUE <- startsWith(colnames(all_data), "LNM_YES_")
index_FALSE <- startsWith(colnames(all_data), "LNM_NO_")
LNM_PRESENT <-subset(all_data, ((startsWith(colnames(all_data), "LNM_YES_")== TRUE)))
LNM_ABSENT <-subset(all_data, ((startsWith(colnames(all_data), "LNM_NO_")== TRUE)))

summary(index_TRUE)
summary(index_FALSE)


LNM_PRESENT$X <- as.list(sapply(LNM_PRESENT$X , function(x) gsub("\'", "", x)))
LNM_ABSENT$X <- as.list(sapply(LNM_ABSENT$X , function(x) gsub("\'", "", x)))
rownames(LNM_PRESENT) <- LNM_PRESENT$X
rownames(LNM_ABSENT) <- LNM_ABSENT$X
LNM_PRESENT$X <- NULL
LNM_ABSENT$X <- NULL

rna_meta_LNM_PRESENT <- LNM_PRESENT[1:5,]
rna_meta_LNM_ABSENT <- LNM_ABSENT[1:5,]

rna_meta_LNM_PRESENT <- as.data.frame(t(rna_meta_LNM_PRESENT))
rna_meta_LNM_ABSENT <- as.data.frame(t(rna_meta_LNM_ABSENT))

rownames(rna_meta_LNM_PRESENT) <- gsub("_", "-", rownames(rna_meta_LNM_PRESENT))
rownames(rna_meta_LNM_ABSENT) <- gsub("_", "-", rownames(rna_meta_LNM_ABSENT))

index <- rna_meta_LNM_PRESENT$`non-cancer cell type` == 0
rna_meta_LNM_PRESENT$`non-cancer cell type`[index] <- "Cancer cell"


index <- rna_meta_LNM_PRESENT$`non-cancer cell type` == "-Fibroblast"
rna_meta_LNM_PRESENT$`non-cancer cell type`[index] <- "Fibroblast"



index <- rna_meta$`non-cancer cell type` == "-Fibroblast"
rna_meta$`non-cancer cell type`[index] <- "Fibroblast"
#Create the counts matrix using the remaining data
rna_clean <- all_data[6:23691,]

names(rna_clean) <- gsub("_", "-", names(rna_clean))
#create seurat objects used for downstream analysis 
rna_seurat <- CreateSeuratObject(rna_clean, project = "puram_data")
rna_seurat <- AddMetaData(rna_seurat, metadata = rna_meta)






all_data$X
all_data$X  <- as.list(sapply(all_data$X , function(x) gsub("\'", "", x)))
rownames(all_data) <- all_data$X
all_data$X <- NULL
#Create rna_meta data using the first 5 rows of the dataframe

rna_meta <- all_data[1:5,]
rna_meta <- as.data.frame(t(rna_meta))
rownames(rna_meta) <- gsub("_", "-", rownames(rna_meta))
index <- rna_meta$`non-cancer cell type` == "-Fibroblast"
rna_meta$`non-cancer cell type`[index] <- "Fibroblast"



#Create the counts matrix using the remaining data
rna_clean <- all_data[6:23691,]

names(rna_clean) <- gsub("_", "-", names(rna_clean))
#create seurat objects used for downstream analysis 
rna_seurat <- CreateSeuratObject(rna_clean, project = "puram_data")
rna_seurat <- AddMetaData(rna_seurat, metadata = rna_meta)


#Perform basic data validation - examine basic characteristics of data
counts_per_cell <- Matrix::colSums(rna_seurat)
counts_per_gene <- Matrix::rowSums(rna_seurat)
genes_per_cell <- Matrix::colSums(rna_seurat)
cells_per_gene <- Matrix::rowSums(rna_seurat)
#Graph the counts per cell, counts per gene, cells per gene, etc 
hist(log10(counts_per_cell+1),main='counts per cell',col='wheat')
hist(log10(genes_per_cell+1), main='genes per cell', col='wheat')
plot(counts_per_cell, genes_per_cell, log='xy', col='wheat')
title('counts vs genes per cell')
hist(log10(counts_per_gene+1), main='counts per gene', col='wheat')
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')

#Normalize the RNA data in the seurat object:
rna_seurat <- NormalizeData(rna_seurat)
rna_seurat <- FindVariableFeatures(rna_seurat)


top10 <- head(VariableFeatures(rna_seurat), 10)
top10

plot1 <- VariableFeaturePlot(rna_seurat)
LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)

rna_seurat <- ScaleData(rna_seurat)

rna_seurat <- RunPCA(rna_seurat, verbose = FALSE)
DimHeatmap(rna_seurat, dims = 1:6, nfeatures = 20, cells = 500, balanced = T)
ElbowPlot(rna_seurat)

rna_seurat <- FindNeighbors(rna_seurat, dims = 1:10)

rna_seurat <- FindClusters(rna_seurat, resolution = 0.8, verbose = FALSE)

rna_seurat <- RunUMAP(rna_seurat, dims = 1:10, verbose = FALSE)
table(rna_seurat@meta.data$seurat_clusters)
DimPlot(rna_seurat, label.size = 4, repel = TRUE, label = TRUE)

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

rna_seurat <- CellCycleScoring(rna_seurat, s.features = s.genes, g2m.features = g2m.genes)
table(rna_seurat[[]]$Phase)

#FeaturePlot(rna_seurat, features = c("LILRA4", "TPM2", "PPBP", "GP1BB"))

rna_seurat <- RunUMAP(rna_seurat, dims = 1:10)
DimPlot(rna_seurat, reduction = "umap", group.by = "non.cancer.cell.type", label = TRUE, repel =TRUE) 

rna_seurat <- RunTSNE(rna_seurat, dims = 1:10)
DimPlot(rna_seurat, reduction = "tsne", group.by = "non.cancer.cell.type", label = TRUE, repel =TRUE)

FeatureScatter(object = rna_seurat, feature1 = "RNF126", feature2 = "PC_1")

VlnPlot(object = rna_seurat, features = c("MIER2"))
