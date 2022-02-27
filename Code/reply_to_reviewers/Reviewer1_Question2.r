library(Seurat)
library(patchwork)
library(ggplot2)
library(harmony)
library(ggpubr)
library(dplyr)
A <- readRDS("Z:/Cailab/Qian_writting/pyroptosis_covid/Integrate/Data/Integrated.rds")
#BALF
BALF <- readRDS('Z:/Cailab/Qian_writting/pyroptosis_covid/Integrate/Data/BALF.rds')
#PBMC
blish <- readRDS('Z:/Cailab/Qian_writting/pyroptosis_covid/Integrate/Data/blish2.rds')
blish$celltype <- blish$cell.type

#CTRL
healthy_lung <- readRDS("Z:/Cailab/Qian_writting/pyroptosis_covid/Integrate/Data/healthy_lung2.rds")
healthy_lung$celltype <- Idents(healthy_lung)
DE_CTRL <- FindAllMarkers(healthy_lung,test.use = 'MAST')
DE_CTRL <- DE_CTRL[DE_CTRL$p_val_adj<0.05,]
DE_CTRL %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5
S2.CTRL <- DotPlot(object = healthy_lung,features = top5$gene, cols = c("lightgrey","red"), group.by = 'celltype')+ RotatedAxis()
write.csv(DE_CTRL,'E:/Pyroptosis_COVID19/reply_to_reviewer/DE_CTRL.csv')
write.csv(top10,'E:/Pyroptosis_COVID19/reply_to_reviewer/DE_CTRL_top10.csv')

PBMC <- subset(A,cells = colnames(A)[colnames(A) %in% colnames(blish)])
Idents(PBMC) <- PBMC$celltype
DE_PBMC <- FindAllMarkers(PBMC,test.use = 'MAST')
DE_PBMC <- DE_PBMC[DE_PBMC$p_val_adj<0.05,]
DE_PBMC %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5
write.csv(DE_PBMC,'E:/Pyroptosis_COVID19/reply_to_reviewer/Reviewer1Q2_DE_PBMC.csv')
write.csv(top10,'E:/Pyroptosis_COVID19/reply_to_reviewer/Reviewer1Q2_DE_PBMC_top10.csv')
S2.PBMC <- DotPlot(object = PBMC,features = top5$gene, cols = c("lightgrey","red"), group.by = 'celltype')+ RotatedAxis()

BALF <- subset(A,cells = colnames(A)[colnames(A) %in% colnames(BALF)])
Idents(BALF) <- BALF$celltype
DE_BALF <- FindAllMarkers(BALF,test.use = 'MAST')
DE_BALF <- DE_BALF[DE_BALF$p_val_adj<0.05,]
DE_BALF <- DE_BALF[!grepl(pattern = 'Antibody-',x = DE_BALF$gene),]
DE_BALF %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5
write.csv(DE_BALF2,'E:/Pyroptosis_COVID19/reply_to_reviewer/Reviewer1Q2_DE_BALF.csv')
write.csv(top10,'E:/Pyroptosis_COVID19/reply_to_reviewer/Reviewer1Q2_DE_BALF_top10.csv')
S2.BALF <- DotPlot(object = BALF,features = top5$gene,cols = c("lightgrey","red"), group.by = 'celltype')+ RotatedAxis()
S2.BALF

png('E:/Pyroptosis_COVID19/reply_to_reviewer/Reviewer1Q2_FigS2.png',res = 800,height = 9000,width = 9000)
S2.CTRL/S2.BALF/S2.PBMC
dev.off()
