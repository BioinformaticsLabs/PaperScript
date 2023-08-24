## Cell-cell interaction analysis

library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(uwot)

sce <- readRDS('sce.rds')

CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB


# Each tissue type for each cancerType
## example: Tumor tissue of UCEC

sce_tumor <- sce[,which(sce$cancerType == "UCEC" & sce$loc == "Tumor")]
cellchat <- createCellChat(object = sce_tumor, group.by = "cluster_name")

cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
write.table(df.net,"UCEC_tumor_net_lr.txt",sep = "\t",row.names = F,quote = F)
cellchat <- computeCommunProbPathway(cellchat)
df.netp <- subsetCommunication(cellchat,slot.name = "netP")

write.table(df.netp,"UCEC_tumor_net_pathway.txt",sep = "\t",row.names = F,quote = F)

cellchat <- aggregateNet(cellchat)

saveRDS(cellchat,file = "UCEC_tumor_cellchat.rds")

# Comparsion among multiple tissue types for each cancer
## example: UCEC

cancer_T <- readRDS("UCEC_tumor_cellchat.rds")
cancer_T <- netAnalysis_computeCentrality(cancer_T,slot.name = "netP")
cancer_T <- computeNetSimilarity(cancer_T, type = "functional")
cancer_T <- netEmbedding(cancer_T, type = "functional", umap.method = "uwot")
cancer_T <- netClustering(cancer_T, type = "functional")

cancer_N <- readRDS("UCEC_normal_cellchat.rds")
cancer_N <- netAnalysis_computeCentrality(cancer_N,slot.name = "netP")
cancer_N <- computeNetSimilarity(cancer_N, type = "functional")
cancer_N <- netEmbedding(cancer_N, type = "functional", umap.method = "uwot")
cancer_N <- netClustering(cancer_N, type = "functional")

cancer_list <- list(Tumor = cancer_T, Normal = cancer_N)

cancer_cellchat <- mergeCellChat(cancer_list,add.names = names(cancer_list),cell.prefix = T)


gg1 <- compareInteractions(cancer_cellchat, show.legend = F, group = c(1,2))+
  scale_fill_manual(values = c("#AD2A47","#6AAD6A"))+
  #scale_fill_manual(values = c("#AD2A47","#009ACD"))+
  #scale_fill_manual(values = c("#AD2A47","#6AAD6A","#009ACD"))+
  theme(axis.title.y=element_text(size=14))+
  theme(axis.text.x = element_text(size = 12,colour = "black"), 
        axis.text.y = element_text(size = 12,colour = "black"))


gg2 <- compareInteractions(cancer_cellchat, show.legend = F, group = c(1,2), measure = "weight")+
  scale_fill_manual(values = c("#AD2A47","#6AAD6A"))+
  #scale_fill_manual(values = c("#AD2A47","#009ACD"))+
  #scale_fill_manual(values = c("#AD2A47","#6AAD6A","#009ACD"))+
  theme(axis.title.y=element_text(size=14))+
  theme(axis.text.x = element_text(size = 12,colour = "black"), 
        axis.text.y = element_text(size = 12,colour = "black"))

pdf("UCEC_Tumor_normal_interaction_num_strength.pdf",width = 5,height = 4)
gg1 + gg2
dev.off()

netVisual_diffInteraction(cancer_cellchat, weight.scale = T, measure = "weight", comparison = c(2,1))

gg1 <- netVisual_heatmap(cancer_cellchat,comparison = c(2,1))
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cancer_cellchat, measure = "weight",comparison = c(2,1))
#> Do heatmap based on a merged object

pdf("UCEC_Tumorvsnormal_Differential_interaction_num_strength.pdf",width = 13,height = 9)
gg1 + gg2
dev.off()

weight.max <- getMaxWeight(cancer_list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(cancer_list)) {
  netVisual_circle(cancer_list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], 
                   edge.width.max = 12, title.name = paste0("Number of interactions - ", names(cancer_list)[i]))
}


num.link <- sapply(cancer_list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(cancer_list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(cancer_list[[i]], title = names(cancer_list)[i], weight.MinMax = weight.MinMax)
}

pdf("UCEC_Tumor_normal_in_out_strength.pdf",width = 13,height = 9)
patchwork::wrap_plots(plots = gg)
dev.off()
