## Data integration and clustering
## gene_count : gene count matrix
## meta : meta data of cells

library(Seurat)
library(harmony)

gene_count <- read.table('/path/gene_count.txt')
meta <- read.table('/path/meta.txt')

sce <- CreateSeuratObject(gene_count,meta.data = meta)

sce <- NormalizeData(sce)

sce <- FindVariableFeatures(sce) 
sce <- ScaleData(sce) %>% RunPCA(verbose=FALSE)

system.time({sce <- RunHarmony(sce, group.by.vars = c("patient_loc"), lambda =2 ,
                                      max.iter.harmony=10,max.iter.cluster = 20)})


sce <- FindNeighbors(sce, reduction = "harmony", dims = 1:30)
sce <- FindClusters(sce, resolution = 1.5)
sce <- RunUMAP(sce, reduction = "harmony", dims = 1:30)
sce <- RunTSNE(sce, reduction = "harmony", dims = 1:30, do.fast=T)

saveRDS(sce,file = "sce.rds")


## Marker genes of each seurat clusters

sce_marker <- FindAllMarkers(sce)
write.table(sce,file = "sce_marker_gene.txt", sep = "\t",row.names = F,quote = F)
