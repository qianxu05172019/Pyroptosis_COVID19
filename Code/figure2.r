A <- readRDS("Z:/Cailab/Qian_writting/pyroptosis_covid/Integrate/Data/Integrated.rds")
##################################################################################################
library(Seurat)
library(patchwork)
library(ggplot2)
library(harmony)
library(ggpubr)
library(dplyr)
library(fgsea)
library(viridis)
library(Nebulosa)
BIOP <- gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=BioPlanet_2019')
table(A$harmony_condition)
table(A$celltype)

#macrophage
macrophage <- A[,A$celltype %in% c('Macrophage')]
macrophage <- NormalizeData(macrophage)
macrophage <- FindVariableFeatures(macrophage)
macrophage <- ScaleData(macrophage)
macrophage <- RunPCA(macrophage)
macrophage <- RunUMAP(macrophage,dims = 1:50)
################ P2A ######################################
P2a <- DimPlot(macrophage,group.by = 'harmony_condition',cols = c("#1b9e77","#7570b3", "#d95f02"))+ggtitle('Condition')+
  xlab(' ')+ylab(' ')+theme(panel.grid = element_blank())+
  theme(axis.text = element_blank())+
  theme(axis.ticks = element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_blank())
P2a
################ P2E ######################################
macrophage <- FindNeighbors(macrophage)
macrophage <- FindClusters(macrophage,resolution = 0.1)
DimPlot(macrophage,label = T,repel = T)+
  xlab(' ')+ylab(' ')+NoLegend()+theme(panel.grid = element_blank())+
  theme(axis.text = element_blank())+
  theme(axis.ticks = element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_blank())
P2E <- DimPlot(macrophage,label = T,repel = T)+
  xlab(' ')+ylab(' ')+theme(panel.grid = element_blank())+
  theme(axis.text = element_blank())+
  theme(axis.ticks = element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_blank())
P2E

################ cluster 3  ######################################
DE <- FindMarkers(macrophage,ident.1 = '3')
DE <- DE[DE$p_val_adj<0.05,]
DE <- DE[!(grepl(pattern = 'Antibody',x = rownames(DE))),]
DE$gene <- rownames(DE)
U <- DE$avg_log2FC
names(U) <- DE$gene
U <- U[!grepl('Rpl|Rps|Rp[[:digit:]]+|Mt-', names(U))]
eD <- fgseaMultilevel(BIOP, U)
eD$leadingEdge <- unlist(lapply(eD$leadingEdge, function(X){paste0(X,collapse = ';')}))
eD <- eD[eD$padj < 0.05,,drop=FALSE]
#write.csv(eD, file = paste0('Table/Macrophage_cluster3_gsea_BIOP.csv')) 

################ P2F ######################################
gSet <- 'NF-kappaB signaling pathway'
pTitle <- 'NF-kappaB signaling\npathway'
P2F <- plotEnrichment(BIOP[[gSet]], U) + 
  labs(title = pTitle) +
  theme_bw() + 
  ylab('Enrichment score') +
  xlab('Gene rank') +
  theme(plot.title = element_text(face = 2, size = 15))
P2F



################ P2G ######################################
gSet <- 'Signal transduction through IL-1R'
pTitle <- 'Signal transduction \nthrough IL-1R'
P2G <- plotEnrichment(BIOP[[gSet]], U) +
  labs(title = pTitle) +
  theme_bw() +
  ylab('Enrichment score') +
  xlab('Gene rank') +
  theme(plot.title = element_text(face = 2, size = 15))
P2G
###############cluster 1#################################

DE <- FindMarkers(macrophage,ident.1 = '1')
DE <- DE[DE$p_val_adj<0.05,]
DE <- DE[!(grepl(pattern = 'Antibody',x = rownames(DE))),]
DE$gene <- rownames(DE)
DE <- DE[DE$p_val_adj<0.05,]
U <- DE$avg_log2FC
names(U) <- DE$gene
U <- U[!grepl('Rpl|Rps|Rp[[:digit:]]+|Mt-', names(U))]
eD <- fgseaMultilevel(BIOP, U)
eD$leadingEdge <- unlist(lapply(eD$leadingEdge, function(X){paste0(X,collapse = ';')}))
eD <- eD[eD$padj < 0.05,,drop=FALSE]
#write.csv(eD, file = paste0('Table/Macrophage_cluster1_gsea_BIOP.csv')) 

gSet <- 'Interferon signaling'
pTitle <- 'Interferon signaling'
#png('Interferon.png',res = 300,height = 1000,width = 1000)
plotEnrichment(BIOP[[gSet]], U) + 
  labs(title = pTitle) +
  theme_bw() + 
  ylab('Enrichment score') +
  xlab('Gene rank') +
  theme(plot.title = element_text(face = 2, size = 15))
#dev.off()

################ P2F ######################################
gSet <- 'Interferon alpha/beta signaling'
pTitle <- 'Interferon alpha/beta \nsignaling'
P2H <- plotEnrichment(BIOP[[gSet]], U) + 
  labs(title = pTitle) +
  theme_bw() + 
  ylab('Enrichment score') +
  xlab('Gene rank') +
  theme(plot.title = element_text(face = 2, size = 15))
P2H

################ P2B ######################################
P2B <- plot_density(macrophage, features = c( 'CASP1', 'IL1B', 'NINJ1'), 
                    joint = TRUE,reduction = 'umap',combine = F)[[4]]+NoLegend()
################ P2C ######################################
P2C <- plot_density(macrophage, features = c( 'MARCO', 'FABP4'), 
                    joint = TRUE,reduction = 'umap',combine = F)[[3]]+NoLegend()

################ P2D ######################################
P2D <- plot_density(macrophage, features = c( 'CXCL10', 'STAT1'), 
                    joint = TRUE,reduction = 'umap',combine = F)[[3]]+NoLegend()


############################ heatmap #############################################
Idents(macrophage) <- macrophage$harmony_condition
macrophage <- FindNeighbors(macrophage)
macrophage <- FindClusters(macrophage,resolution = 0.1)
DE <- FindAllMarkers(macrophage)
DE <- DE[DE$p_val_adj <0.05&DE$avg_log2FC >2,]
DE <- DE[!(grepl(pattern = 'Antibody|ATP|MALAT1|MT|RPS|RPL',x = rownames(DE))),]
top <- DE %>% group_by(cluster) %>% top_n(50, avg_log2FC)

P2I <- DotPlot(object = macrophage, features  = top$gene,cols = c('#377eb8','#e41a1c'))+
  theme(axis.text.x = element_text(angle =90,size = 10))
P2I 

P2B <- P2B+ggtitle('Severe: CASP1+IL1B+NINJ1')
P2D <- P2D+ ggtitle('Moderate: CXCL10+STAT1')
P2C <- P2C + ggtitle('Healthy: MARCO+FABP4')
P2E <- P2E+NoLegend()

png('Z:/Cailab/Qian_writting/pyroptosis_covid/Integrate/Figure/Fig2.png',res = 500,height = 6000,width = 9000)
(P2a+(P2B|P2D|P2C) + 
    plot_layout(widths = c(1, 4))) /
  (P2E+ P2I + 
  plot_layout(widths = c(1, 4)))/
  (P2F|P2G|P2H)
dev.off()
