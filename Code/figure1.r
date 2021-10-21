setwd("Z:/Cailab/Qian_writting/pyroptosis_covid/Integrate")
library(Seurat)
library(patchwork)
library(ggplot2)
library(harmony)
library(ggpubr)
library(dplyr)
A <- readRDS("Z:/Cailab/Qian_writting/pyroptosis_covid/Integrate/Integrated.PBMC.BALF.HealthyLung.rds")
DimPlot(A)


A$Experiment
BALF <- subset(A,cells = colnames(A)[A$Experiment %in% c('BALF')])
fig1a1 <- DimPlot(BALF,group.by = 'harmony_cell_type')+xlab('UMAP 1')+ylab('UMAP 2')+ggtitle('')+ theme_classic()
fig1a1

DimPlot(A,group.by = 'harmony_cell_type',label = T,repel = T)
A <- FindNeighbors(A)
A <- FindClusters(A,resolution = 0.05)
DimPlot(A,label = T,repel = T)
all.markers <- FindAllMarkers(object = A)
Top20 <- all.markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)
write.csv(Top20,'DE_top20.all.cluster.csv')

plots <- VlnPlot(A, features = c("MS4A1",'CD79A',#b cell
                                 'CD2','CD3G', # CD4 T
                                 'IL1R2',# cDC
                                 'CD24',#epithelial
                                 'AQP9','LGMN',#macrophage
                                 'OAS1','CD86','CCR2',# monocyte
                                 'FUT4',#neutrophil
                                 'GNLY'), #NK
                group.by = 'harmony_cell_type',
                 pt.size = 0, combine = FALSE)
png('Z:/Cailab/Qian_writting/pyroptosis_covid/Integrate/celltype_markers.png',res = 400,height = 7000,width = 7000)
wrap_plots(plots = plots, ncol = 2)
dev.off()
FeaturePlot(A,c('FUT4'))


B <- A
DimPlot(B)
B <- NormalizeData(B)
B <- FindVariableFeatures(B)
B <- ScaleData(B)
B <- RunPCA(B)
B <- RunUMAP(B,dims = 1:20)
DimPlot(B)
DimPlot(B,group.by = 'singler')
Dimplot(B,group.by = 'singler')
B <- FindNeighbors(B)
B <- FindClusters(B,resolution = 0.05)
DimPlot(B)
all.markers <- FindAllMarkers(object = B)
Top20 <- all.markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)
write.csv(Top20,'DE_top20.clustersofmacrophagelike.csv')
table(Idents(B))


setwd("Z:/Cailab/Qian_writting/pyroptosis_covid/Integrate")
A <- readRDS("Integrated.PBMC.BALF.HealthyLung.rds")
library(SingleR)
library(celldex)
library(SingleCellExperiment)
ref <- BlueprintEncodeData()
sce <- A@assays$RNA@data
pretend.cell.labels <- sample(letters, ncol(sce), replace=TRUE)
pretend.gene.lengths <- sample(10000, nrow(sce),replace=TRUE)

sce <- SingleCellExperiment(list(counts=A@assays$RNA@counts,logcounts = A@assays$RNA@data),
                            colData=DataFrame(label=pretend.cell.labels),
                            rowData=DataFrame(length=pretend.gene.lengths),
                            metadata=list(study="Integrated")
)

pred <- SingleR(test=sce, ref=ref, labels=ref$label.main)
table(pred$labels)

A$singler <- pred$labels[match(colnames(A),rownames(pred))]
DimPlot(A,group.by = 'singler',label = T)
table(A$singler)
A <- A[,A$singler %in% c('DC','Macrophages','Monocytes','Neutrophils','B-cells','CD4+ T-cells',
                         'CD8+ T-cells','Epithelial cells','NK cells')]
DimPlot(A,group.by = 'singler',label = T)
DimPlot(A,group.by = 'Experiment')

plots <- VlnPlot(A, features = c("CASP1",'GSDMD','IL18',
                                 'IL1B','CASP8',
                                 'CASP4','NINJ1','TLR4','NLRP3',
                                 'NLRC4','ACTB','PPIA')
                 ,split.by = "harmony_condition",group.by = 'harmony_cell_type',
                 pt.size = 0, combine = FALSE)
png('Z:/Cailab/Qian_writting/pyroptosis_covid/Integrate/singleR/py_markers_harmony_Ventilated_1.png',res = 400,height = 7000,width = 7000)
wrap_plots(plots = plots, ncol = 3)
dev.off()
pred$scores

