library(Seurat)
# Dataset downloaded from file:///Z:/Cailab/Qian_writting/pyroptosis_covid/reply_to_reviewer/suplementaryData1/lambrechtslab%20-%20immune%20atlas.html
# Step 1 prepare dataset
Allcells.counts <- readRDS("Z:/Cailab/Qian_writting/pyroptosis_covid/reply_to_reviewer/suplementaryData1/2069-Allcells.counts.rds")
meta <- read.csv('Z:/Cailab/Qian_writting/pyroptosis_covid/reply_to_reviewer/suplementaryData1/2080-Allcells.meta.data.csv')

colnames(Allcells.counts)
table(meta$CellType)
meta <- meta[meta$CellType %in% c('Alveolar_macrophage','B_cell','CD4_Tcell','CD8_Tcell','cDC','Md_macrophage','Monocyte','NK','pDC'),]
Allcells.counts <- Allcells.counts[,colnames(Allcells.counts) %in% meta$Cell]

All <- CreateSeuratObject(Allcells.counts)
All$CellType <- meta$CellType[match(colnames(All),meta$Cell)]
table(All$CellType)
All$Desease <- meta$Disease[match(colnames(All),meta$Cell)]
All$PatientType <- meta$PatientType[match(colnames(All),meta$Cell)]

All <- NormalizeData(All)
All <- FindVariableFeatures(All)
All <- ScaleData(All)
All <- RunPCA(All)
All <- RunUMAP(All,dims=1:50)
DimPlot(All,group.by = 'PatientType',label = T, repel = T)
All$CellType
DimPlot(All,group.by = 'Desease',label = T, repel = T)
saveRDS(All, 'Z:/Cailab/Qian_writting/pyroptosis_covid/reply_to_reviewer/suplementaryData1/Reviewer2.dataset1.rds')
