## ConsensusClusterPlus and survival analysis

library(ConsensusClusterPlus)
library(pheatmap)
library(survival)
library(survminer)

gene_id <- read.table("gencode.v35lift37.annotation.txt",sep = "\t",header = T,stringsAsFactors = F)
gene_id_new <- unique(gene_id[,c(1,7)])

regulon14 <- gene_id_new[which(gene_id_new$Hugo_Symbol %in% c("ETV1","FOXP3","ETV7","RUNX2","IRF4","PRDM1","EZH2",
                                                              "FLI1","ELF2","IKZF1","KLF2","KLF3","KLF13","TBX21")),]

rownames(regulon14) <- regulon14$Gene_ID


tcga_cancer <- c('ACC','BLCA','BRCA','CESC','CHOL','COAD','DLBC','ESCA','GBM','HNSC','KICH','KIRC','KIRP','LAML','LGG','LIHC',
                 'LUAD','LUSC','MESO','OV','PAAD','PCPG','PRAD','READ','SARC','SKCM','STAD','TGCT','THCA','THYM','UCEC','UCS','UVM')


for(i in tcga_cancer){
  
  dir.create(paste0("path",i))
  
  setwd(paste0("path",i))
  
  tumor_sample <- read.table(paste0("path",i,"_tumor.txt"),sep = "\t",stringsAsFactors = F)

  tcga_cli <- read.table(paste0("path",i,".txt"),sep = "\t",header = T,row.names = 1,stringsAsFactors = F)
  
  tcga_cli_tumor <- tcga_cli[unique(tumor_sample$V1),]
  
  tcga_exp <- read.table(paste0("path",i,".txt"),sep = "\t",header = T,row.names = 1,stringsAsFactors = F)
  
  colnames(tcga_exp) <- gsub("\\.","-",colnames(tcga_exp))
  
  
  tumor_exp <- tcga_exp[intersect(regulon14$Gene_ID,rownames(tcga_exp)),tumor_sample$V1]
  
  rownames(tumor_exp) <- regulon14[rownames(tumor_exp),"Hugo_Symbol"]
  tumor_exp <- as.matrix(tumor_exp)

  
  tumor_exp_scale <- t(scale(t(tumor_exp)))
  
  P <- ConsensusClusterPlus(tumor_exp_scale,maxK=3,reps=1000,pItem = 0.8,
                            pFeature = 0.8,
                            distance = "spearman",
                            clusterAlg="pam",
                            plot="pdf",
                            writeTable=T)
  
  
  
  exp_quantile <- tumor_exp
  #exp_quantile <- as.matrix(exp_quantile)
  
  
  for(m in rownames(exp_quantile)){
    
    tmp_exp <- exp_quantile[m,]
    
    tmp_quantile <- quantile(tmp_exp,c(0.25,0.5,0.75))
    
    
    tmp25 <- tmp_exp[tmp_exp < tmp_quantile[1]]
    tmp5 <- tmp_exp[tmp_exp >= tmp_quantile[1] & tmp_exp < tmp_quantile[2]]
    tmp75 <- tmp_exp[tmp_exp >= tmp_quantile[2] & tmp_exp < tmp_quantile[3]]
    tmp1 <- tmp_exp[tmp_exp >= tmp_quantile[3]]
    
    
    exp_quantile[m,names(tmp25)] <- 0.25
    exp_quantile[m,names(tmp5)] <- 0.5
    exp_quantile[m,names(tmp75)] <- 0.75
    exp_quantile[m,names(tmp1)] <- 1
    
    
  }
  
  sample_class <- read.csv("untitled_consensus_cluster/untitled_consensus_cluster.k=2.consensusClass.csv",header = F,stringsAsFactors = F)
  
  sample_class <- sample_class[order(sample_class$V2),]
  sample_class$V2 <- paste0("class_",sample_class$V2)
  
  
  colnames(sample_class) <- c("sample","class")
  rownames(sample_class) <- sample_class$sample
  
  
  exp_quantile_new <- exp_quantile[,sample_class$sample]
  
  
  annotation_col <- data.frame(group = factor(sample_class$class))
  ann_colors <- list(group = c(class_1= "#8ED4C8",class_2= "#FFFFB5"))
  rownames(annotation_col) <- colnames(exp_quantile_new)
  
  P_exp<- pheatmap(exp_quantile_new,
                   cluster_cols = F, cluster_rows = T,
                   annotation_col= annotation_col, 
                   annotation_colors= ann_colors,
                   show_colnames = F, 
                   border_color = NA,
                   col=colorRampPalette(c("#93C6DF","#D2E6F0","#FDDCC8","#F4A783"))(4))
  
  
  pdf(paste0(i,"_exp_level_pheatmap.pdf"),width = 5,height = 3)
  print(P_exp)
  dev.off()
  
  
  tcga_cli_tumor$class <- sample_class[rownames(tcga_cli_tumor),"class"]
  
  fit_mean <- survfit(Surv(OS.time,OS)~class,data=tcga_cli_tumor)
  
  P_sur <- ggsurvplot(fit_mean,risk.table="abs_pct",#生存统计统计衿
                      
                      conf.int=F,#添加置信区间帿
                      
                      palette = c("#F39F6F","#679C6F"),#颜色设置
                      
                      pval=TRUE,#log-rank检骿
                      
                      pval.method=TRUE)+
    labs(x="Survival time (days)")
  
  pdf(paste0(i,"_sig_gene_survivalplot.pdf"),width = 5,height = 5)
  print(P_sur,newpage = F)
  dev.off()


}
