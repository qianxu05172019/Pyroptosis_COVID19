setwd("Z:/Cailab/Qian_writting/pyroptosis_covid/PBMCs")
library(patchwork)
library(Seurat)
blish_covid.seu <- readRDS("Z:/Cailab/Qian_writting/pyroptosis_covid/PBMCs/blish_covid.seu.rds")
table(blish_covid.seu$singler)
table(blish_covid.seu$cell.type)

blish_bcell <-  blish_covid.seu[,blish_covid.seu$cell.type %in% c('B')]
set.seed(111)
bcell <- sample(x=colnames(blish_bcell),size = 2000,replace = F)
blish_bcell  <- blish_bcell[,colnames(blish_bcell) %in% bcell]
table(blish_bcell$Ventilated)

blish_NK <- blish_covid.seu[,blish_covid.seu$cell.type %in% c('NK')]
set.seed(111)
cell <- sample(x=colnames(blish_NK),size = 2000,replace = F)
blish_NK  <- blish_NK[,colnames(blish_NK) %in% cell]
table(blish_NK$Ventilated)

blish_Monocyte <- blish_covid.seu[,blish_covid.seu$cell.type %in% c('CD14 Monocyte')]
set.seed(111)
cell <- sample(x=colnames(blish_Monocyte),size = 2000,replace = F)
blish_Monocyte  <- blish_Monocyte[,colnames(blish_Monocyte) %in% cell]
#blish_Monocyte$cell.type <- 'Monocyte'
table(blish_Monocyte$Ventilated)
blish_dc <- blish_covid.seu[,blish_covid.seu$cell.type %in% c('DC','pDC')]
blish_dc$cell.type <- 'DC'

blish_other <- blish_covid.seu[,blish_covid.seu$cell.type %in% c('CD4 T','Neutrophil')]
table(blish_other$Ventilated)


blish <- merge(blish_bcell,blish_NK)
blish <- merge(blish,blish_Monocyte)
blish <- merge(blish,blish_other)
blish <- merge(blish,blish_dc)

table(blish$singler)
#saveRDS(blish,file = 'Z:/Cailab/Qian_writting/pyroptosis_covid/Integrate/Data/blish.rds')
#blish <- readRDS('Z:/Cailab/Qian_writting/pyroptosis_covid/Integrate/Data/blish.rds')
table(blish$singler)
table(blish$Status)
table(blish$Ventilated)


plots <- VlnPlot(blish, features = c("CASP1",'GSDMD','IL18',
                                 'IL1B',
                                 'CASP4','NINJ1','TLR4',
                                 'PYCARD',
                                 'GAPDH','PPIA')
                 ,group.by = 'cell.type',split.by = 'Ventilated',
                 pt.size = 0, combine = FALSE)
png('Z:/Cailab/Qian_writting/pyroptosis_covid/PBMCs/py_markertest.png',res = 400,height = 7000,width = 7000)
wrap_plots(plots = plots, ncol =3)
dev.off()

blish$cell.type[blish$cell.type %in% c('B')] <- 'B cell'
blish$cell.type[blish$cell.type %in% c('CD14 Monocyte')] <- 'Monocyte'
blish$cell.type[blish$cell.type %in% c('CD4 T')] <- 'CD4+ T-cell'
blish$cell.type[blish$cell.type %in% c('DC')] <- 'cDC'
blish$cell.type[blish$cell.type %in% c('Neutrophil')] <- 'Neutrophil'
blish$cell.type[blish$cell.type %in% c('NK')] <- 'NK cell'


saveRDS(blish,file = 'Z:/Cailab/Qian_writting/pyroptosis_covid/Integrate/Data/blish2.rds')
