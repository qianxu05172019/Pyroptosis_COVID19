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

# FeaturePlot(macrophage,features = c('NINJ1','CASP1','IL1B',
#                                     'CASP4','MARCO','FABP4',
#                                     'IL1R1','NLRP3','GSDMD'))

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
  theme(axis.text.x = element_text(angle =45))
P2I 

png('Figure/Fig2.png',res = 600,height = 6000,width = 9000)
(P2a|P2B|P2C|P2D) /
  (P2E|P2F|P2G|P2H)/
  P2I
dev.off()


########################### Generate correlation of each pair of genes #############################
df <- A[rownames(A) %in% c('CASP1', 'PYCARD','IL1B','GSDMD','NINJ1','ACTB','PPIA','GAPDH','LDHA'),]
df <- as.data.frame(df@assays$RNA@data)
df <- t(df)
df <- as.data.frame(df)
df$Group <- A$harmony_condition
Group <- as.factor(df$Group)
df$color <- NA
df$color[df$Group %in% c('severe')] <- '#d95f02'
df$color[df$Group %in% c('moderate')] <- '#7570b3'
df$color[df$Group %in% c('Healthy')] <- '#1b9e77'

df$size <- NA
df$size[df$Group %in% c('severe')] <- 1
df$size[df$Group %in% c('moderate')] <- 0.8
df$size[df$Group %in% c('Healthy')] <- 0.5

df$alpha <- NA
df$alpha[df$Group %in% c('severe')] <- 1
df$alpha[df$Group %in% c('moderate')] <- 0.5
df$alpha[df$Group %in% c('Healthy')] <- 0.4


df <- df[df$CASP1>0,]
df <- df[df$LDHA>0,]
png('Cor_CASP1_LDHA.png',res = 300,height = 1000,width = 1000)
ggplot(data = df, aes(x =CASP1, 
                      y = LDHA)) +
  geom_point(alpha = df$alpha,colour = df$color,size=df$size)+
  scale_color_discrete("Group") +theme_bw()+stat_cor(method="pearson")
dev.off()

plots <- VlnPlot(A, features = c('NLRP1', 'NLRP2', 'NLRP3','NLRC4',  'AIM2','GSDMD','LDHA','GAPDH','CASP3','IL18'), 
                 group.by ='celltype', split.by = 'harmony_condition',
                 pt.size = 0, combine = FALSE)
png('Z:/Cailab/Qian_writting/pyroptosis_covid/Integrate/test.png',res = 400,height = 7000,width = 7000)
wrap_plots(plots = plots, ncol = 3)
dev.off()

