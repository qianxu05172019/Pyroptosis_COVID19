library(Seurat)
library(patchwork)
library(ggplot2)
library(readxl)
library(dplyr)
setwd("Z:/Cailab/Qian_writting/pyroptosis_covid/")
A <- readRDS('reply_to_reviewer/Data/BALF2.rds')
meta <- read.csv('reply_to_reviewer/suplementaryData1/2080-Allcells.meta.data.patientID_added.csv')
rownames(meta) <- meta$Cell
table(A$Condition)
n <- A@assays$RNA@data
n <- n[rownames(n) %in% c("TLR1",'TLR2','TLR3',
                          'TLR4','TLR5','TLR6',
                          'TLR7','TLR8','TLR9',
                          'TLR10',
                          'NLRP1','NLRP2','NLRP3','NLRP4',
                          'NLRP5','NLRP6','NLRP7','NLRP8',
                          'NLRP9','NLRP10','NLRP11','NLRP12',
                          'NLRP13','NLRP14'),]
n <- as.data.frame(t(as.matrix(n)))
n$Condition <- A$Condition
n$celltype <- A$CellType
n$cellType_condition <- paste0(n$celltype,'_',n$Condition)
table(n$cellType_condition)
n <- n[!n$Condition %in% 'Control',]

library(ggpubr)
CellType <- n$celltype
Condition <- n$Condition
TLR1 <- ggviolin(n, x = "celltype", y = "TLR1", color = "Condition", 
                      add = "jitter", legend = "none",size=0,fill = 'Condition',palette = c("#7570b3", "#d95f02"),xlab = '')+ 
  stat_compare_means(aes(group = Condition), label = "p.signif",hide.ns = T)+
  rotate_x_text(angle = 45)
TLR2 <- ggviolin(n, x = "celltype", y = "TLR2", color = "Condition", 
                 add = "jitter", legend = "none",size=0,fill = 'Condition',palette = c("#7570b3", "#d95f02"),xlab = '')+ 
  stat_compare_means(aes(group = Condition), label = "p.signif",hide.ns = T)+
  rotate_x_text(angle = 45)
TLR3 <- ggviolin(n, x = "celltype", y = "TLR3", color = "Condition", 
                 add = "jitter", legend = "none",size=0,fill = 'Condition',palette = c("#7570b3", "#d95f02"),xlab = '')+ 
  stat_compare_means(aes(group = Condition), label = "p.signif",hide.ns = T)+
  rotate_x_text(angle = 45)
TLR4 <- ggviolin(n, x = "celltype", y = "TLR4", color = "Condition", 
                 add = "jitter", legend = "none",size=0,fill = 'Condition',palette = c("#7570b3", "#d95f02"),xlab = '')+ 
  stat_compare_means(aes(group = Condition), label = "p.signif",hide.ns = T)+
  rotate_x_text(angle = 45)
TLR5 <- ggviolin(n, x = "celltype", y = "TLR5", color = "Condition", 
                 add = "jitter", legend = "none",size=0,fill = 'Condition',palette = c("#7570b3", "#d95f02"),xlab = '')+ 
  stat_compare_means(aes(group = Condition), label = "p.signif",hide.ns = T)+
  rotate_x_text(angle = 45)
TLR6 <- ggviolin(n, x = "celltype", y = "TLR6", color = "Condition", 
                 add = "jitter", legend = "none",size=0,fill = 'Condition',palette = c("#7570b3", "#d95f02"),xlab = '')+ 
  stat_compare_means(aes(group = Condition), label = "p.signif",hide.ns = T)+
  rotate_x_text(angle = 45)
TLR7 <- ggviolin(n, x = "celltype", y = "TLR7", color = "Condition", 
                 add = "jitter", legend = "none",size=0,fill = 'Condition',palette = c("#7570b3", "#d95f02"),xlab = '')+ 
  stat_compare_means(aes(group = Condition), label = "p.signif",hide.ns = T)+
  rotate_x_text(angle = 45)
TLR8 <- ggviolin(n, x = "celltype", y = "TLR8", color = "Condition", 
                 add = "jitter", legend = "none",size=0,fill = 'Condition',palette = c("#7570b3", "#d95f02"),xlab = '')+ 
  stat_compare_means(aes(group = Condition), label = "p.signif",hide.ns = T)+
  rotate_x_text(angle = 45)
TLR9 <- ggviolin(n, x = "celltype", y = "TLR9", color = "Condition", 
                 add = "jitter", legend = "none",size=0,fill = 'Condition',palette = c("#7570b3", "#d95f02"),xlab = '')+ 
  stat_compare_means(aes(group = Condition), label = "p.signif",hide.ns = T)+
  rotate_x_text(angle = 45)
TLR10 <- ggviolin(n, x = "celltype", y = "TLR10", color = "Condition", 
                 add = "jitter", legend = "none",size=0,fill = 'Condition',palette = c("#7570b3", "#d95f02"),xlab = '')+ 
  stat_compare_means(aes(group = Condition), label = "p.signif",hide.ns = T)+
  rotate_x_text(angle = 45)


NLRP1 <- ggviolin(n, x = "celltype", y = "NLRP1", color = "Condition", 
                  add = "jitter", legend = "none",size=0,fill = 'Condition',palette = c("#7570b3", "#d95f02"),xlab = '')+ 
  stat_compare_means(aes(group = Condition), label = "p.signif",hide.ns = T)+
  rotate_x_text(angle = 45)
NLRP2 <- ggviolin(n, x = "celltype", y = "NLRP2", color = "Condition", 
                  add = "jitter", legend = "none",size=0,fill = 'Condition',palette = c("#7570b3", "#d95f02"),xlab = '')+ 
  stat_compare_means(aes(group = Condition), label = "p.signif",hide.ns = T)+
  rotate_x_text(angle = 45)
NLRP3 <- ggviolin(n, x = "celltype", y = "NLRP3", color = "Condition", 
                  add = "jitter", legend = "none",size=0,fill = 'Condition',palette = c("#7570b3", "#d95f02"),xlab = '')+ 
  stat_compare_means(aes(group = Condition), label = "p.signif",hide.ns = T)+
  rotate_x_text(angle = 45)
NLRP4 <- ggviolin(n, x = "celltype", y = "NLRP4", color = "Condition", 
                  add = "jitter", legend = "none",size=0,fill = 'Condition',palette = c("#7570b3", "#d95f02"),xlab = '')+ 
  stat_compare_means(aes(group = Condition), label = "p.signif",hide.ns = T)+
  rotate_x_text(angle = 45)
NLRP5 <- ggviolin(n, x = "celltype", y = "NLRP5", color = "Condition", 
                  add = "jitter", legend = "none",size=0,fill = 'Condition',palette = c("#7570b3", "#d95f02"),xlab = '')+ 
  stat_compare_means(aes(group = Condition), label = "p.signif",hide.ns = T)+
  rotate_x_text(angle = 45)
NLRP6 <- ggviolin(n, x = "celltype", y = "NLRP6", color = "Condition", 
                  add = "jitter", legend = "none",size=0,fill = 'Condition',palette = c("#7570b3", "#d95f02"),xlab = '')+ 
  stat_compare_means(aes(group = Condition), label = "p.signif",hide.ns = T)+
  rotate_x_text(angle = 45)
NLRP7 <- ggviolin(n, x = "celltype", y = "NLRP7", color = "Condition", 
                  add = "jitter", legend = "none",size=0,fill = 'Condition',palette = c("#7570b3", "#d95f02"),xlab = '')+ 
  stat_compare_means(aes(group = Condition), label = "p.signif",hide.ns = T)+
  rotate_x_text(angle = 45)
NLRP8 <- ggviolin(n, x = "celltype", y = "NLRP8", color = "Condition", 
                  add = "jitter", legend = "none",size=0,fill = 'Condition',palette = c("#7570b3", "#d95f02"),xlab = '')+ 
  stat_compare_means(aes(group = Condition), label = "p.signif",hide.ns = T)+
  rotate_x_text(angle = 45)
NLRP9 <- ggviolin(n, x = "celltype", y = "NLRP9", color = "Condition", 
                  add = "jitter", legend = "none",size=0,fill = 'Condition',palette = c("#7570b3", "#d95f02"),xlab = '')+ 
  stat_compare_means(aes(group = Condition), label = "p.signif",hide.ns = T)+
  rotate_x_text(angle = 45)
NLRP10 <- ggviolin(n, x = "celltype", y = "NLRP10", color = "Condition", 
                   add = "jitter", legend = "none",size=0,fill = 'Condition',palette = c("#7570b3", "#d95f02"),xlab = '')+ 
  stat_compare_means(aes(group = Condition), label = "p.signif",hide.ns = T)+
  rotate_x_text(angle = 45)
NLRP11 <- ggviolin(n, x = "celltype", y = "NLRP10", color = "Condition", 
                   add = "jitter", legend = "none",size=0,fill = 'Condition',palette = c("#7570b3", "#d95f02"),xlab = '')+ 
  stat_compare_means(aes(group = Condition), label = "p.signif",hide.ns = T)+
  rotate_x_text(angle = 45)
NLRP12 <- ggviolin(n, x = "celltype", y = "NLRP10", color = "Condition", 
                   add = "jitter", legend = "none",size=0,fill = 'Condition',palette = c("#7570b3", "#d95f02"),xlab = '')+ 
  stat_compare_means(aes(group = Condition), label = "p.signif",hide.ns = T)+
  rotate_x_text(angle = 45)
NLRP13 <- ggviolin(n, x = "celltype", y = "NLRP10", color = "Condition", 
                   add = "jitter", legend = "none",size=0,fill = 'Condition',palette = c("#7570b3", "#d95f02"),xlab = '')+ 
  stat_compare_means(aes(group = Condition), label = "p.signif",hide.ns = T)+
  rotate_x_text(angle = 45)
NLRP14 <- ggviolin(n, x = "celltype", y = "NLRP10", color = "Condition", 
                   add = "jitter", legend = "none",size=0,fill = 'Condition',palette = c("#7570b3", "#d95f02"),xlab = '')+ 
  stat_compare_means(aes(group = Condition), label = "p.signif",hide.ns = T)+
  rotate_x_text(angle = 45)

png('Z:/Cailab/Qian_writting/pyroptosis_covid/reply_to_reviewer/Figure/Supplementary Fig. S1.png',res = 600,height = 9000,width = 9000)
(TLR1|TLR2|TLR3|TLR4|TLR5)  /
  (TLR6|TLR7|TLR8|TLR9|TLR10) /
  (NLRP1|NLRP2|NLRP3|NLRP4|NLRP5)  /
  (NLRP6|NLRP7|NLRP8|NLRP9|NLRP10) /
  (NLRP11|NLRP12|NLRP13|NLRP14)
dev.off()
