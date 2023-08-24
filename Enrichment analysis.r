## Enrichment analysis

library(clusterprofiler)
library(simplifyEnrichment)
library(GOSemSim)
library(org.Hs.eg.db)


CD8_marker <- read.table('CD8_marker.txt',sep='\t',header=T,stringsAsFactors=F)
TXNIP_marker <- CD8_marker[which(CD8_marker$cluster_name == "CD8_TXNIP+ T" & CD8_marker$avg_logFC > 0),"gene"]


data(geneList, package="DOSE")
TXNIP_marker_gene = bitr(TXNIP_marker,
                fromType="SYMBOL", 
                toType="ENTREZID",  
                OrgDb="org.Hs.eg.db") 

ego_ALL <- enrichGO(gene = TXNIP_marker_gene$ENTREZID, 
                    universe = names(geneList), 
                    OrgDb = org.Hs.eg.db, 
                    #keytype = 'ENSEMBL',
                    ont = "ALL", 
                    pAdjustMethod = "BH", 
                    pvalueCutoff = 1, 
                    qvalueCutoff = 1,
                    readable = TRUE) 


write.table(as.data.frame(ego_ALL),
            "TXNIP_markerGene_clusterprofiler_goenrich.txt",
            sep = "\t",quote = F)


kegg_ALL <- enrichKEGG(gene = TXNIP_marker_gene$ENTREZID, 
                       organism="hsa",
                       keyType = 'kegg',
                       #ont = "ALL", 
                       pAdjustMethod = "BH", 
                       pvalueCutoff = 1, 
                       qvalueCutoff = 1) 


write.table(as.data.frame(kegg_ALL),
            "TXNIP_markerGene_clusterprofiler_KEGGenrich.txt",
            sep = "\t",quote = F)



TXNIP_GO <- read.table('TXNIP_markerGene_clusterprofiler_goenrich.txt',
                       sep = '\t',header = T,row.names = 1,stringsAsFactors = F,quote = '')

TXNIP_GO <- TXNIP_GO[which(TXNIP_GO$ONTOLOGY == 'BP' & TXNIP_GO$p.adjust < 0.05),]

matBP <- GO_similarity(TXNIP_GO$ID,ont='BP',db=org.Hs.eg.db)

df <- simplifyGO(matBP,plot=T,fontsize_range=c(12,18),method='walktrap')
