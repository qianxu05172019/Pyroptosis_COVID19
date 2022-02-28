library(Seurat)
library(patchwork)
library(ggplot2)
library(harmony)
library(ggpubr)
library(dplyr)
library(rstatix)
setwd("Z:/Cailab/Qian_writting/pyroptosis_covid/reply_to_reviewer")
BALF <- readRDS('Z:/Cailab/Qian_writting/pyroptosis_covid/Integrate/Data/BALF.rds')
table(BALF$orig.ident)
table(BALF$condition)
table(BALF$celltype)


healthy_lung <- readRDS("Z:/Cailab/Qian_writting/pyroptosis_covid/Integrate/Data/healthy_lung2.rds")
healthy_lung$celltype <- Idents(healthy_lung)

blish <- readRDS('Z:/Cailab/Qian_writting/pyroptosis_covid/Integrate/Data/blish2.rds')
table(blish$Ventilated)
blish$celltype <- blish$cell.type
table(blish$celltype)

BALF2 <- readRDS('Z:/Cailab/Qian_writting/pyroptosis_covid/reply_to_reviewer/Data/BALF2.rds')
meta <- read.csv('Z:/Cailab/Qian_writting/pyroptosis_covid/reply_to_reviewer/suplementaryData1/2080-Allcells.meta.data.patientID_added.csv')
table(meta$PatientNumber)
rownames(meta) <- meta$X
table(meta$CellType)
BALF2 <- BALF2[,colnames(BALF2) %in% meta$Cell[meta$CellType %in% c('B_cell', 'CD4_Tcell','CD8_Tcell','cDC',
                                                                    'Md_macrophage','Monocyte','Neutrophil','NK')]]
dim(BALF2)
table(BALF2$CellType)

################################# load data #############################################################
A <- merge(healthy_lung,BALF)
A <- merge(A,blish)
A <- merge(A,BALF2)
dim(A)

A <- NormalizeData(A)
A <- FindVariableFeatures(A)
A <- ScaleData(A)
A <- RunPCA(A)
A <- RunUMAP(A,dims = 1:20)

table(A$orig.ident)
################# set the batch ###################################################################
A$Experiment <- NA
A$Experiment[A$orig.ident %in% c('COV002','COV004','COV007','COV012','COV034','COV036')] <- 'BALF1'
A$Experiment[A$orig.ident %in% c("covid_555_1","covid_555_2","covid_556",
                                 "covid_557" , "covid_558" , "covid_559" ,
                                 "covid_560" , "covid_561",  "HIP002" ,
                                 "HIP015" , "HIP023","HIP043",
                                 "HIP044","HIP045")] <- 'PBMC'
A$Experiment[A$orig.ident %in% c("Macrophage","NK cell", "CD4+ T cell" , "CD8+ T cell","B cell",
                                 "Monocyte" ,"cDC")] <- 'Healthy lung'
A$Experiment[A$orig.ident %in% unique(meta$PatientNumber)] <- 'BALF2'

############################# set the cell type ##################################################
table(A$celltype)
dim(A)
A$celltype[A$celltype %in% c('NK','NK_cell','NK cell')] <- 'NK'
A$celltype[A$celltype %in% c('Neutrophils','Neutrophil')] <- 'Neutrophil'
A$celltype[A$celltype %in% c('CD4+ T-cell','CD4+ T cell')] <- 'CD4+ T cell'
A$celltype[A$celltype %in% c('CD8+ T-cell','CD8+ T cell')] <- 'CD8+ T cell'
A$celltype[A$celltype %in% c('B cell','B_cell')] <- 'B cell'

################ BALF2 dataset ###############
A$celltype[A$CellType %in% c('B_cell')] <- 'B cell'
A$celltype[A$CellType %in% c('CD4_Tcell')] <- 'CD4+ T cell'
A$celltype[A$CellType %in% c('CD8_Tcell')] <- 'CD8+ T cell'
A$celltype[A$CellType %in% c('cDC')] <- 'cDC'
A$celltype[A$CellType %in% c('Md_macrophage')] <- 'Macrophage'
A$celltype[A$CellType %in% c('Monocyte')] <- 'Monocyte'
A$celltype[A$CellType %in% c('Neutrophil')] <- 'Neutrophil'
A$celltype[A$CellType %in% c('NK')] <- 'NK'

A$harmony_condition <- NA
A$harmony_condition[A$orig.ident %in% c('COV007','COV034','COV036')] <- 'Severe'
A$harmony_condition[A$orig.ident %in% c('COV002','COV004','COV012')] <- 'Moderate'
A$harmony_condition[A$Ventilated %in% c('Healthy')] <- 'Healthy'
A$harmony_condition[A$Ventilated %in% c('NonVent')] <- 'Moderate'
A$harmony_condition[A$Ventilated %in% c('Vent')] <- 'Severe'
A$harmony_condition[A$orig.ident %in% c("Macrophage",
                                        "NK cell","CD4+ T-cell",
                                        "CD8+ T-cell" ,"B cell" ,
                                        "cDC" )] <- 'Healthy'
A$harmony_condition[A$Condition %in% c('Moderate')] <- 'Moderate'
A$harmony_condition[A$Condition %in% c('Severe')] <- 'Severe'

table(A$harmony_condition)
A <- RunHarmony(A,group.by.vars = 'Experiment')
A <- RunUMAP(A,reduction = 'harmony',dims = 1:20)
DimPlot(A,group.by = 'harmony_condition')
saveRDS(A,"Z:/Cailab/Qian_writting/pyroptosis_covid/Integrate/Data/Integrated_4datagroup.rds")

A <- readRDS("Z:/Cailab/Qian_writting/pyroptosis_covid/Integrate/Data/Integrated_4datagroup.rds")
dim(A)
##################################################### P1A1 ######################################
P1A1 <- DimPlot(A,group.by = 'celltype',label = T,repel = T,
                cols = c('#d53e4f','#f46d43','#fdae61','#fee08b','#e6f598',
                         '#abdda4','#66c2a5','#3288bd'))+xlab(' ')+ylab(' ')+
  ggtitle('Integrated')
P1B1 <- DimPlot(A,group.by = 'harmony_condition',cols = c("#1b9e77","#7570b3", "#d95f02"))+xlab(' ')+ylab(' ')+ggtitle('Integrated')
P1B1
################ P1A2 ######################################
PBMC <- subset(A,cells = colnames(A)[colnames(A) %in% colnames(blish)])
P1A2 <- DimPlot(PBMC,group.by = 'celltype',label = T,repel = T,
                cols = c('#d53e4f','#f46d43' ,'#fee08b' ,'#abdda4','#66c2a5','#3288bd'))+xlab(' ')+ylab(' ')+ggtitle('PBMC')+NoLegend()+theme(panel.grid = element_blank())+
  theme(axis.text = element_blank())+
  theme(axis.ticks = element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_blank())

P1B2 <- DimPlot(PBMC,group.by = 'harmony_condition',cols = c("#1b9e77","#7570b3", "#d95f02"))+NoLegend()+xlab(' ')+ylab(' ')+
  ggtitle('PBMC')+NoLegend()+theme(panel.grid = element_blank())+
  theme(axis.text = element_blank())+
  theme(axis.ticks = element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_blank())
P1B2
################ P1A3 ######################################
BALF <- subset(A,cells = colnames(A)[colnames(A) %in% colnames(BALF)])
P1A3 <- DimPlot(BALF,group.by = 'celltype',label = T,repel = T,
                cols = c(        '#f46d43','#fdae61','#fee08b',
                                 '#e6f598',      '#66c2a5','#3288bd'))+
  xlab(' ')+ylab(' ')+ggtitle('BALF1')+NoLegend()+theme(panel.grid = element_blank())+
  theme(axis.text = element_blank())+
  theme(axis.ticks = element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_blank())


P1B3 <- DimPlot(BALF,group.by = 'harmony_condition',cols = c( "#7570b3", "#d95f02"))+
  xlab(' ')+ylab(' ')+ggtitle('BALF1')+theme(panel.grid = element_blank())+NoLegend()+
  theme(axis.text = element_blank())+
  theme(axis.ticks = element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_blank())
P1B3 
############### P1A4 ######################################
healthy_lung <- subset(A,cells = colnames(A)[colnames(A) %in% colnames(healthy_lung)])
P1A4 <- DimPlot(healthy_lung,group.by = 'celltype',label = T,repel = T,
                cols = c('#d53e4f','#f46d43','#fdae61', '#e6f598','#3288bd'))+xlab(' ')+ylab(' ')+NoLegend()+
  ggtitle('Healthy lung')+NoLegend()+theme(panel.grid = element_blank())+
  theme(axis.text = element_blank())+
  theme(axis.ticks = element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_blank())

P1B4 <- DimPlot(healthy_lung,group.by = 'harmony_condition',
                cols = c("#1b9e77"))+xlab(' ')+ylab(' ')+
  ggtitle('Healthy lung')+NoLegend()+theme(panel.grid = element_blank())+
  theme(axis.text = element_blank())+
  theme(axis.ticks = element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_blank())
P1B4

################ P1A5 ######################################
BALF2 <- subset(A,cells = colnames(A)[colnames(A) %in% colnames(BALF2)])
P1A5 <- DimPlot(BALF2,group.by = 'celltype',label = T,repel = T,
                cols = c( c('#d53e4f','#f46d43','#fdae61','#fee08b','#e6f598',
                            '#abdda4','#66c2a5','#3288bd')))+
  xlab(' ')+ylab(' ')+ggtitle('BALF2')+NoLegend()+theme(panel.grid = element_blank())+
  theme(axis.text = element_blank())+
  theme(axis.ticks = element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_blank())


P1B5 <- DimPlot(BALF2,group.by = 'harmony_condition',cols = c( "#7570b3", "#d95f02"))+
  xlab(' ')+ylab(' ')+ggtitle('BALF2')+theme(panel.grid = element_blank())+NoLegend()+
  theme(axis.text = element_blank())+
  theme(axis.ticks = element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_blank())
P1B5 

#############################

A <- readRDS("Z:/Cailab/Qian_writting/pyroptosis_covid/Integrate/Data/Integrated_4datagroup.rds")
table(A$CellType)
n <- A@assays$RNA@data
n <- n[rownames(n) %in% c("CASP1",'GSDMD','IL18',
                          'IL1B','PYCARD',
                          'CASP4','NINJ1','TLR2',
                          'NLRC5'),]
table(A$harmony_condition)
n <- as.data.frame(t(as.matrix(n)))
n$Condition <- A$harmony_condition
n$celltype <- A$celltype
CellType <- n$celltype
Condition <- n$Condition
R2Q1_IL1B <- ggviolin(n, x = "celltype", y = "IL1B", color = "Condition", ylim=c(0,6), 
                      add = "jitter", legend = "none",size=0,add.params = list(size=0.1),fill = 'Condition',palette = c("#1b9e77","#7570b3", "#d95f02"),xlab = '')+
  rotate_x_text(angle = 45)

R2Q1_GSDMD <- ggviolin(n, x = "celltype", y = "GSDMD", color = "Condition", ylim=c(0,6), 
                       add = "jitter", legend = "none",size=0, add.params = list(size=0.1),fill = 'Condition',palette = c("#1b9e77","#7570b3", "#d95f02"),xlab = '')+
  rotate_x_text(angle = 45)

R2Q1_NINJ1 <- ggviolin(n, x = "celltype", y = "NINJ1", color = "Condition", ylim=c(0,6), 
                       add = "jitter", legend = "none",size=0, add.params = list(size=0.1),fill = 'Condition',palette = c("#1b9e77","#7570b3", "#d95f02"),xlab = '')+
  rotate_x_text(angle = 45)

R2Q1_TLR2 <- ggviolin(n, x = "celltype", y = "TLR2", color = "Condition", ylim=c(0,6), 
                      add = "jitter", legend = "none",size=0, add.params = list(size=0.1),fill = 'Condition',palette = c("#1b9e77","#7570b3", "#d95f02"),xlab = '')+
  rotate_x_text(angle = 45)
R2Q1_TLR2

R2Q1_CASP4 <- ggviolin(n, x = "celltype", y = "CASP4", color = "Condition", ylim=c(0,6), 
                       add = "jitter", legend = "none",size=0, add.params = list(size=0.1),fill = 'Condition',palette = c("#1b9e77","#7570b3", "#d95f02"),xlab = '')+
  rotate_x_text(angle = 45)

R2Q1_IL18 <- ggviolin(n, x = "celltype", y = "IL18", color = "Condition", ylim=c(0,6),
                      add = "jitter", legend = "none",size=0, add.params = list(size=0.1),fill = 'Condition',palette = c("#1b9e77","#7570b3", "#d95f02"),xlab = '')+
  rotate_x_text(angle = 45)

R2Q1_CASP1 <- ggviolin(n, x = "celltype", y = "CASP1", color = "Condition", ylim=c(0,6), 
                       add = "jitter", legend = "none",size=0, add.params = list(size=0.1),fill = 'Condition',palette = c("#1b9e77","#7570b3", "#d95f02"),xlab = '')+
  rotate_x_text(angle = 45)

R2Q1_NLRC5 <- ggviolin(n, x = "celltype", y = "NLRC5", color = "Condition", ylim=c(0,6), 
                       add = "jitter", legend = "none",size=0,add.params = list(size=0.1),fill = 'Condition',palette = c("#1b9e77","#7570b3", "#d95f02"),xlab = '')+
  rotate_x_text(angle = 45)
R2Q1_NLRC5

png('Z:/Cailab/Qian_writting/pyroptosis_covid/reply_to_reviewer/Figure/Fig1.png',res = 700,height = 10000,width = 12000)
(P1A1|P1A4|P1A2|P1A3|P1A5)/
  (P1B1|P1B4|P1B2|P1B3|P1B5)/
  (R2Q1_CASP1|R2Q1_CASP4|R2Q1_IL1B|R2Q1_IL18)  /
  (R2Q1_NINJ1|R2Q1_GSDMD|R2Q1_NLRC5|R2Q1_TLR2)
dev.off()
