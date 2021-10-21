library(Seurat)
library(patchwork)
library(ggplot2)
library(harmony)
all <- readRDS("Z:/Cailab/Qian_writting/pyroptosis_covid/BALF/all.rds")
table(all$celltype)
all$condition<- NA
all$condition[all$orig.ident %in% c('COV007','COV034','COV036')] <- 'severe'
all$condition[all$orig.ident %in% c('COV002','COV004','COV012')] <- 'moderate'

all_CD4_moderate <- all[,all$celltype %in% c('CD4+ T-cell') & all$condition %in% c('moderate')]
set.seed(111)
cell <- sample(colnames(all_CD4_moderate),size = 400,replace = F)
all_CD4_moderate <- all_CD4_moderate[,colnames(all_CD4_moderate) %in% cell]
table(all_CD4_moderate$condition)

all_CD4_severe <- all[,all$celltype %in% c('CD4+ T-cell')& all$condition %in% c('severe')]
set.seed(111)
cell <- sample(colnames(all_CD4_severe),size = 400,replace = F)
all_CD4_severe <- all_CD4_severe[,colnames(all_CD4_severe) %in% cell]
table(all_CD4_severe$condition)


all_CD8_moderate <- all[,all$celltype %in% c('CD8+ T-cell') & all$condition %in% c('moderate')]
set.seed(111)
cell <- sample(colnames(all_CD8_moderate),size = 400,replace = F)
all_CD8_moderate <- all_CD8_moderate[,colnames(all_CD8_moderate) %in% cell]
table(all_CD8_moderate$condition)

all_CD8_severe <- all[,all$celltype %in% c('CD8+ T-cell') & all$condition %in% c('severe')]
set.seed(111)
cell <- sample(colnames(all_CD8_severe),size = 400,replace = F)
all_CD8_severe <- all_CD8_severe[,colnames(all_CD8_severe) %in% cell]
table(all_CD8_severe$condition)

all_Neutrophil_severe <- all[,all$celltype %in% c('Neutrophil') & all$orig.ident %in% c('COV007','COV034','COV036')]
dim(all_Neutrophil_severe)
set.seed(111)
cell <- sample(colnames(all_Neutrophil_severe),size = 400,replace = F)
all_Neutrophil_severe <- all_Neutrophil_severe[,colnames(all_Neutrophil_severe) %in% cell]
table(all_Neutrophil_severe$condition)

all_Neutrophil_moderate <- all[,all$celltype %in% c('Neutrophil') & all$orig.ident %in% c('COV002','COV004','COV012')]
dim(all_Neutrophil_moderate)
cell <- sample(colnames(all_Neutrophil_moderate),size = 400,replace = F)
all_Neutrophil_moderate <- all_Neutrophil_moderate[,colnames(all_Neutrophil_moderate) %in% cell]
table(all_Neutrophil_moderate$condition)


all_Macrophage_severe <- all[,all$celltype %in% c('Macrophage') & all$orig.ident %in% c('COV007','COV034','COV036')]
dim(all_Macrophage_severe)
cell <- sample(colnames(all_Macrophage_severe),size = 400,replace = F)
all_Macrophage_severe <- all_Macrophage_severe[,colnames(all_Macrophage_severe) %in% cell]
table(all_Macrophage_severe$condition)



all_Macrophage_moderate <- all[,all$celltype %in% c('Macrophage') & all$orig.ident %in% c('COV002','COV004','COV012')]
dim(all_Macrophage_moderate)
set.seed(111)
cell <- sample(colnames(all_Macrophage_moderate),size = 400,replace = F)
all_Macrophage_moderate <- all_Macrophage_moderate[,colnames(all_Macrophage_moderate) %in% cell]
table(all_Macrophage_moderate$condition)


all_NK_moderate <- all[,all$celltype %in% c('NK') & all$orig.ident %in% c('COV002','COV004','COV012')]
dim(all_NK_moderate)
set.seed(111)
cell <- sample(colnames(all_NK_moderate),size = 400,replace = F)
all_NK_moderate <- all_NK_moderate[,colnames(all_NK_moderate) %in% cell]
table(all_NK_moderate$condition)

all_NK_severe <- all[,all$celltype %in% c('NK') & all$condition %in% c('severe')]
dim(all_NK_severe)

table(all$celltype)
all_other <- all[,all$celltype %in% c('cDC')]

all <- merge(all_CD4_moderate,all_CD4_severe)
all <- merge(all,all_CD8_moderate)
all <- merge(all,all_CD8_severe)
all <- merge(all,all_Neutrophil_severe)
all <- merge(all,all_Neutrophil_moderate)
all <- merge(all,all_Macrophage_severe)
all <- merge(all,all_Macrophage_moderate)
all <- merge(all,all_NK_moderate)
all <- merge(all,all_NK_severe)
all <- merge(all,all_other)
dim(all)
table(all$celltype)
A <- all
A <- NormalizeData(A)
A <- FindVariableFeatures(A)
A <- ScaleData(A)
A <- RunPCA(A)
A <- RunUMAP(A,dims=1:50)
DimPlot(A,label = T,repel = T,group.by = 'celltype')
A <- RunHarmony(A,group.by.vars = 'orig.ident')
A <- RunUMAP(A,reduction = 'harmony',dims = 1:20)
DimPlot(A,label = T,repel = T,group.by = 'celltype')
DimPlot(A,label = T,repel = T,group.by = 'condition')
plots <- VlnPlot(A, features = c("CASP1",'GSDMD','IL18',
                                 'IL1B',
                                 'CASP4','NINJ1','TLR4',
                                 'PYCARD',
                                 'GAPDH','PPIA')
                 ,group.by = 'celltype',split.by = 'condition',
                 pt.size = 0, combine = FALSE)
png('Z:/Cailab/Qian_writting/pyroptosis_covid/BALF/py_marker.png',res = 400,height = 4000,width = 7000)
wrap_plots(plots = plots, ncol = 5)
dev.off()

saveRDS(A,file = 'Z:/Cailab/Qian_writting/pyroptosis_covid/Integrate/Data/BALF.rds')

