library(Seurat)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(matrixStats)
computPyroptosis <- function(sampleSize,SeuratObject,marker) {
  b <- SeuratObject@assays$RNA@counts
  sd <- rowSds(as.matrix(b))
  c <- as.data.frame(cbind(rownames(b),rowMeans(b),sd))
  colnames(c) <- c('gene','rowmeans','SD')
  c$rowmeans <- as.numeric(c$rowmeans)
  c <- c[order(c$rowmeans,decreasing = F),]
  c$logCV2 <- log((as.numeric(c$SD)/(c$rowmeans))^2)
  c$logMean <- log(c$rowmeans)
  TID <- c()
  ID <- c()
  set.seed(1)
  for (i in 1:length(marker)) {
    y <- c$logCV2[c$gene %in% marker[i]]
    x <- c$logMean[c$gene %in% marker[i]]
    TID[i] <- which(c$logCV2==y &c$logMean==x)
    pool <- which(c$logCV2 < y+1 & c$logCV2 >y-1 & c$logMean > x-0.02 & c$logMean < x+0.02)
    ID[i] <- sample(x = pool,size = 1,replace = T)
    message('Sampling from ',length(pool),' genes with similar mean and logCV2')
  }
  gene_star <- rownames(c)[ID]
  c[rownames(c) %in%c(gene_star,marker),]
  delta <- c()
  selected_cell <- sample(colnames(b),size = sampleSize,replace = F)
  subset_cell <- subset(SeuratObject,cells = selected_cell)
  count_subset_cell <- as.data.frame(subset_cell@assays$RNA@counts)
  
  for (k in 1:dim(count_subset_cell)[2]) {
    cell <- as.data.frame(count_subset_cell[,k])
    colnames(cell) <- colnames(count_subset_cell)[k]
    score_star <- sum(cell[which(rownames(count_subset_cell) %in% gene_star),])
    score <- sum(cell[which(rownames(count_subset_cell) %in% marker),])
    delta[k] <- score-score_star
    message('Processing cell ', k, ' of ', dim(count_subset_cell)[2])
  }
  f <- as.data.frame(cbind(colnames(count_subset_cell),delta))
  subset_cell$delta <- f$delta[match(colnames(subset_cell),f$V1)]
  Score <- as.numeric(subset_cell$delta)
  Condition <- subset_cell$harmony_condition
  CellType <- subset_cell$celltype
  table(subset_cell$celltype)
  subset_cell$Phagocytes <- NA 
  subset_cell$Phagocytes[which(subset_cell$celltype %in% c('Macrophage','Monocyte','Neutrophil','cDC'))] <- 'Phagocytes'
  subset_cell$Phagocytes[which(!subset_cell$celltype %in% c('Macrophage','Monocyte','Neutrophil','cDC'))] <- 'Lymphocytes'
  Phagocytes <- subset_cell$Phagocytes
  df <- as.data.frame(cbind(Condition,Score,CellType,Phagocytes))
  df$Score <- as.numeric(df$Score)
  
  f_Condition <- factor(Condition,levels = c('Healthy','moderate','severe'))
  df$f_Condition <- f_Condition
  df <- df[order(df$Score,decreasing = T),]
  return(df)
}
A <- readRDS("Z:/Cailab/Qian_writting/pyroptosis_covid/Integrate/Data/Integrated.rds")
marker <- c("CASP1",'IL1B','NINJ1','PYCARD','CASP4','TLR2','GSDMD','NLRC5')
df <- computPyroptosis(sampleSize = 5000,SeuratObject = A,marker = marker)

ggboxplot(df, "Phagocytes", "Score",
                 fill = "Condition", palette = c("#1b9e77","#7570b3", "#d95f02"),
                 outlier.shape = NA)+ylim(-10,10)+xlab('')+theme_bw()+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

