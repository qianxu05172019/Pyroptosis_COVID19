library(anndata)
library(Seurat)
library(scTenifoldNet)
library(patchwork)
library(harmony)
library(ggplot2)
setwd("C:/Users/qianxu0517/Downloads")
library(reticulate)
use_python("C:/Users/qianxu0517/software/anaconda3/python.exe")
ad <- read_h5ad("BALF_VIB-UGent_processed_cleaned.h5ad")
ad <- ad$X
dim(ad)
rownames(ad)[1:10]

# Patient COV002, COV004, and COV012, who was discharged alive were considered to be in moderate condition
rm(a1)
a1 <- ad[grep(pattern = 'COV002_',rownames(ad)),]
a1 <- a1[,which(colSums(a1)>50)]
dim(a1)
saveRDS(a1,'Z:/Cailab/Qian_writting/pyroptosis_covid/BALF/COV002.rds')

rm(a1)
a1 <- ad[grep(pattern = 'COV004_',rownames(ad)),]
a1 <- a1[,which(colSums(a1)>50)]
dim(a1)
saveRDS(a1,'Z:/Cailab/Qian_writting/pyroptosis_covid/BALF/COV004.rds')

rm(a1)
a1 <- ad[grep(pattern = 'COV012_',rownames(ad)),]
a1 <- a1[,which(colSums(a1)>50)]
dim(a1)
saveRDS(a1,'Z:/Cailab/Qian_writting/pyroptosis_covid/BALF/COV012.rds')

# Patient COV007, COV034, and COV036, who was dead from COVID19 were considered to be in severe condition
rm(a1)
a1 <- ad[grep(pattern = 'COV007_',rownames(ad)),]
a1 <- a1[,which(colSums(a1)>50)]
dim(a1)
saveRDS(a1,'Z:/Cailab/Qian_writting/pyroptosis_covid/BALF/Death_COV007.rds')

rm(a1)
a1 <- ad[grep(pattern = 'COV034_',rownames(ad)),]
a1 <- a1[,which(colSums(a1)>50)]
dim(a1)
saveRDS(a1,'Z:/Cailab/Qian_writting/pyroptosis_covid/BALF/Death_COV034.rds')

rm(a1)
a1 <- ad[grep(pattern = 'COV036_',rownames(ad)),]
a1 <- a1[,which(colSums(a1)>50)]
dim(a1)
saveRDS(a1,'Z:/Cailab/Qian_writting/pyroptosis_covid/BALF/Death_COV036.rds')

ad <- read_h5ad("BALF_VIB-UGent_processed_cleaned.h5ad")
metadata <- ad$obs
write.csv(metadata,file = 'Z:/Cailab/Qian_writting/pyroptosis_covid/BALF/metadata.csv')

setwd("Z:/Cailab/Qian_writting/pyroptosis_covid/BALF")

meta <- read.csv('metadata.csv')
COV002 <- readRDS("Z:/Cailab/Qian_writting/pyroptosis_covid/BALF/COV002.rds")
COV002 <- t(as.matrix(COV002))
COV002 <- CreateSeuratObject(scQC(COV002),project ='alive_COV002')
COV002$celltype <- meta$Annotation[match(colnames(COV002),meta$X)]

COV012 <- readRDS("Z:/Cailab/Qian_writting/pyroptosis_covid/BALF/COV012.rds")
COV012 <- t(as.matrix(COV012))
COV012 <- CreateSeuratObject(scQC(COV012),project ='alive_COV012')
COV012$celltype <- meta$Annotation[match(colnames(COV012),meta$X)]

COV004 <- readRDS("Z:/Cailab/Qian_writting/pyroptosis_covid/BALF/COV004.rds")
COV004 <- t(as.matrix(COV004))
COV004 <- CreateSeuratObject(scQC(COV004),project ='alive_COV004')
COV004$celltype <- meta$Annotation[match(colnames(COV004),meta$X)]

A <- merge(COV002,COV004)
rm(COV002,COV004)
A <- merge(A,COV012)
rm(COV012)

Death_COV034 <- readRDS("Z:/Cailab/Qian_writting/pyroptosis_covid/BALF/Death_COV034.rds")
Death_COV034 <- t(as.matrix(Death_COV034))
Death_COV034 <- CreateSeuratObject(scQC(Death_COV034,minLibSize = 500),project ='Death_COV034')
Death_COV034$celltype <- meta$Annotation[match(colnames(Death_COV034),meta$X)]
A <- merge(A,Death_COV034)
rm(Death_COV034)

Death_COV036 <- readRDS("Z:/Cailab/Qian_writting/pyroptosis_covid/BALF/Death_COV036.rds")
Death_COV036 <- t(as.matrix(Death_COV036))
Death_COV036 <- CreateSeuratObject(scQC(Death_COV036,minLibSize = 500),project ='Death_COV036')
Death_COV036$celltype <- meta$Annotation[match(colnames(Death_COV036),meta$X)]
A <- merge(A,Death_COV036)
rm(Death_COV036)

Death_COV007 <- readRDS("Z:/Cailab/Qian_writting/pyroptosis_covid/BALF/Death_COV007.rds")
Death_COV007 <- t(as.matrix(Death_COV007))
Death_COV007 <- CreateSeuratObject(scQC(Death_COV007,minLibSize = 500),project ='Death_COV007')
Death_COV007$celltype <- meta$Annotation[match(colnames(Death_COV007),meta$X)]
A <- merge(A,Death_COV007)
rm(Death_COV007)

A[["percent.mt"]] <- PercentageFeatureSet(A, pattern = "^MT-")
VlnPlot(A, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

A <- NormalizeData(A)
A <- FindVariableFeatures(A)
A <- ScaleData(A)
A <- RunPCA(A)
A <-  FindNeighbors(A)
A <- FindClusters(A,resolution = 0.1)
A <- RunUMAP(A,dims = 1:20)
saveRDS(A,file='Z:/Cailab/Qian_writting/pyroptosis_covid/BALF/all.rds')
A <- readRDS('Z:/Cailab/Qian_writting/pyroptosis_covid/BALF/all.rds')
UMAPPlot(A,group.by='orig.ident')+xlab('UMAP 1')+ylab('UMAP 2')+theme_bw()


A <- RunHarmony(A,group.by.vars = 'orig.ident')
A <- RunUMAP(A,reduction = 'harmony',dims = 1:20)

DimPlot(A)
DimPlot(A,group.by = 'orig.ident')
A$condition <- NA
A$condition[which(A$orig.ident %in% c('COV002','COV004','COV012'))] <- 'Recovered'
A$condition[which(A$orig.ident %in% c('COV007','COV034','COV036'))] <- 'Dead'
DimPlot(A,group.by = 'condition')+xlab('UMAP 1')+ylab('UMAP 2')+theme_bw()
