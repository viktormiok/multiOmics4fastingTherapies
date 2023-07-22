rm(list=ls())

# load libraries
pkgs = c("DESeq2",
         "ggplot2",
         "gplots",
         "RColorBrewer",
         "clusterProfiler",
         "org.Mm.eg.db",
         "stringdist",
         "dplyr",
         "ggrepel") # package names
inst <- suppressMessages(lapply(pkgs, 
                                library,
                                character.only=TRUE)
)

load("/Users/viktorian.miok/Documents/consultation/Katarina/Integration_Data/RNA_seq_data/Initial_Analysis.RData")

sampleName <- read.table(file='~/Documents/consultation/Katarina/Integration_Data/RNA_seq_data/Sample_names.txt',
                         header=T,
                         sep='\t'
)
colnames(expression.matrix) = as.character(sampleName[,2])
options(repr.plot.width=10, repr.plot.height=10)

# Obsolete terms
noliv <- read.csv("/Users/viktorian.miok/Documents/consultation/Katarina/rna_protein_integration/input/200828 non-liver related terms.csv",
                  header=FALSE
)
obsolete <- read.csv("/Users/viktorian.miok/Documents/consultation/Katarina/rna_protein_integration/input/200728_obsolete_GO_terms.csv",
                     header=FALSE
)

###############################################################################
#
#   TRF study
#
###############################################################################

trfdat = expression.matrix[,c(64:71,56:58,60:63)]

id = colnames(trfdat)
condition = as.factor(c(rep("TRFHFD",8), rep("ALHFD",7)))
metaData = data.frame(id, condition)

dds <- DESeqDataSetFromMatrix(countData=trfdat,
                              colData=metaData,
                              design=~condition
)
featureData=data.frame(gene=rownames(dds))
mcols(dds)=DataFrame(mcols(dds), featureData)

design(dds) =~ condition  
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds,
                           fit='parametric'
)

keep = rowSums(counts(dds)) >= 10
dds = dds[keep,]

dds <- DESeq(dds)
res_alhfd_vs_trfhfd = results(dds)

res = res_alhfd_vs_trfhfd[order(res_alhfd_vs_trfhfd$pvalue),]
resultsTRF = as.data.frame(dplyr::mutate(as.data.frame(res),
                                         significant=ifelse(res$pvalue < 0.05, "pvalue<0.05", "not signif.")),
                           row.names=rownames(res))
resultsTRF = resultsTRF[(!is.na(resultsTRF$pvalue)),]

p <- ggplot2::ggplot(resultsTRF, ggplot2::aes(log2FoldChange, -log10(pvalue))) +
     ggplot2::geom_point(ggplot2::aes(col=significant)) + 
     ggplot2::scale_color_manual(values=c("grey", "dodgerblue3")) +
     guides(col=guide_legend(nrow=2)) + 
     ggplot2::ggtitle("") + 
     xlim(-5,5) + 
     ylim(0,20)
p

vsd <- assay(vst(dds, 
                 fitType="local",
                 blind=TRUE)
)

df_pca = prcomp(t(vsd))
df_out = as.data.frame(df_pca$x)

p <- ggplot(df_out,aes(x=PC1,y=PC2,color=condition)) +
     ggplot2::scale_color_manual(values=c("dodgerblue4", "dodgerblue1")) +
     geom_point(size=5)
p

#######################################
# UP
#######################################
sig.gene <- bitr(rownames(resultsTRF[which(resultsTRF$sig == "FDR<0.05" & resultsTRF$log2FoldChange > 0),]),
                 fromType="ENSEMBL",
                 toType="ENTREZID",
                 OrgDb=org.Mm.eg.db
)
kk1 <- enrichKEGG(gene=sig.gene[,2],
                  organism='mmu'
)
ego1 <- enrichGO(gene=sig.gene[,2], 
                 OrgDb=org.Mm.eg.db, 
                 ont="BP", 
                 pAdjustMethod="BH",
                 readable=TRUE
)
ego2 <- enrichGO(gene=sig.gene[,2],
                 OrgDb=org.Mm.eg.db, 
                 ont="CC",
                 pAdjustMethod="BH",
                 readable=TRUE
)
ego3 <- enrichGO(gene=sig.gene[,2], 
                 OrgDb=org.Mm.eg.db,
                 ont="MF",
                 pAdjustMethod="BH", 
                 readable=TRUE
)

bp1 = ego1@result
setwd('/Users/viktorian.miok/Documents/consultation/Katarina/Integration_Data/RNA_seq_data/results/TRF')
#write.csv(bp1,"TRF_GO_BP_up_RNA.csv")

bp1 = bp1[bp1$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp1[,4]))
bp1 = bp1[(lk > 15) & (lk < 500),]
bp1 = bp1[bp1$Count > 4,]
bp1 = bp1[!bp1$ID %in% noliv[,1],]
bp1 = bp1[!bp1$ID %in% obsolete[,1],]
bp1 = bp1[order(bp1$p.adjust,decreasing=TRUE),]
#write.csv(bp1,"TRF_GO_BP_trimed_up_RNA_a.csv")

rk = character()
for(i in 1:(nrow(bp1)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp1)){
    pk = c(pk,sum(unlist(strsplit(bp1[i,8], "/")) %in% unlist(strsplit(bp1[j,8], "/")))/length(unlist(strsplit(bp1[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk,"ostavi")
bp1 = bp1[rk1 == "ostavi",]
#write.csv(bp1,"TRF_GO_BP_trimed_up_RNA_b.csv")

bp2 = ego2@result
#write.csv(bp2,"TRF_GO_CC_up_RNA.csv")

bp2 = bp2[bp2$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp2[,4]))
bp2 = bp2[(lk > 15) & (lk < 500),]
bp2 = bp2[bp2$Count > 4,]
bp2 = bp2[!bp2$ID %in% noliv[,1],]
bp2 = bp2[!bp2$ID %in% obsolete[,1],]
bp2 = bp2[order(bp2$p.adjust,decreasing=TRUE),]
#write.csv(bp2,"TRF_GO_CC_trimed_up_RNA_a.csv")

rk = character()
for(i in 1:(nrow(bp2)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp2)){
    pk = c(pk,sum(unlist(strsplit(bp2[i,8], "/")) %in% unlist(strsplit(bp2[j,8], "/")))/length(unlist(strsplit(bp2[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk,"ostavi")
bp2 = bp2[rk1 == "ostavi",]
#write.csv(bp2,"TRF_GO_CC_trimed_up_RNA_b.csv")

bp3 = ego3@result
#write.csv(bp3,"TRF_GO_MF_up_RNA.csv")

bp3 = bp3[bp3$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp3[,4]))
bp3 = bp3[(lk > 15) & (lk < 500),]
bp3 = bp3[bp3$Count > 4,]
bp3 = bp3[!bp3$ID %in% noliv[,1],]
bp3 = bp3[!bp3$ID %in% obsolete[,1],]
bp3 = bp3[order(bp3$p.adjust,decreasing=TRUE),]
#write.csv(bp3,"TRF_GO_MF_trimed_up_RNA_a.csv")

rk = character()
for(i in 1:(nrow(bp3)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp3)){
    pk = c(pk,sum(unlist(strsplit(bp3[i,8], "/")) %in% unlist(strsplit(bp3[j,8], "/")))/length(unlist(strsplit(bp3[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk,"ostavi")
bp3 = bp3[rk1 == "ostavi",]
#write.csv(bp3,"TRF_GO_MF_trimed_up_RNA_b.csv")

bp4 = kk1@result
for(i in 1:nrow(bp4)){
  bp4[i,8] <- paste(suppressMessages(bitr(unlist(strsplit(bp4[i,8], "/")),
                                          fromType="ENTREZID",
                                          toType="SYMBOL",
                                          OrgDb=org.Mm.eg.db)[,2]),
                    collapse="/"
)
}

#write.csv(bp4,"TRF_KEGG_up_RNA.csv")

bp4 = bp4[bp4$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp4[,4]))
bp4 = bp4[(lk > 15) & (lk < 500),]
bp4 = bp4[bp4$Count > 4,]
bp4 = bp4[!bp4$ID %in% noliv[,1],]
bp4 = bp4[!bp4$ID %in% obsolete[,1],]
bp4 = bp4[order(bp4$p.adjust,decreasing=TRUE),]

#write.csv(bp4,"TRF_KEGG_trimed_up_RNA_a.csv")

rk = character()
for(i in 1:(nrow(bp4)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp4)){
    pk = c(pk,sum(unlist(strsplit(bp4[i,8], "/")) %in% unlist(strsplit(bp4[j,8], "/")))/length(unlist(strsplit(bp4[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk,"ostavi")
bp4 = bp4[rk1 == "ostavi",]
#write.csv(bp4,"TRF_KEGG_trimed_up_RNA_b.csv")

#######################################
# DOWN
#######################################

sig.gene <- bitr(rownames(resultsTRF[which(resultsTRF$sig == "FDR<0.05" & resultsTRF$log2FoldChange < 0),]),
                 fromType="ENSEMBL",
                 toType="ENTREZID",
                 OrgDb=org.Mm.eg.db
)
kk1 <- enrichKEGG(gene=sig.gene[,2],
                  organism='mmu'
)
ego1 <- enrichGO(gene=sig.gene[,2], 
                 OrgDb=org.Mm.eg.db, 
                 ont="BP",
                 pAdjustMethod="BH", 
                 readable=TRUE
)
ego2 <- enrichGO(gene=sig.gene[,2], 
                 OrgDb=org.Mm.eg.db,
                 ont="CC", 
                 pAdjustMethod="BH",
                 readable=TRUE
)
ego3 <- enrichGO(gene=sig.gene[,2],
                 OrgDb=org.Mm.eg.db,
                 ont="MF",
                 pAdjustMethod="BH",
                 readable=TRUE
)

bp1 = ego1@result
setwd('/Users/viktorian.miok/Documents/consultation/Katarina/Integration_Data/RNA_seq_data/results/TRF')
#write.csv(bp1,"TRF_GO_BP_down_RNA.csv")

bp1 = bp1[bp1$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp1[,4]))
bp1 = bp1[(lk > 15) & (lk < 500),]
bp1 = bp1[bp1$Count > 4,]
bp1 = bp1[!bp1$ID %in% noliv[,1],]
bp1 = bp1[!bp1$ID %in% obsolete[,1],]
bp1 = bp1[order(bp1$p.adjust,decreasing=TRUE),]# (as.numeric(sub("\\/.*", "", bp1[,4]))
#write.csv(bp1,"TRF_GO_BP_trimed_down_RNA_a.csv")

rk = character()
for(i in 1:(nrow(bp1)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp1)){
    pk = c(pk,sum(unlist(strsplit(bp1[i,8], "/")) %in% unlist(strsplit(bp1[j,8], "/")))/length(unlist(strsplit(bp1[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk,"ostavi")
bp1 = bp1[rk1 == "ostavi",]
#write.csv(bp1,"TRF_GO_BP_trimed_down_RNA_b.csv")

bp2 = ego2@result
#write.csv(bp2,"TRF_GO_CC_down_RNA.csv")

bp2 = bp2[bp2$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp2[,4]))
bp2 = bp2[(lk > 15) & (lk < 500),]
bp2 = bp2[bp2$Count > 4,]
bp2 = bp2[!bp2$ID %in% noliv[,1],]
bp2 = bp2[!bp2$ID %in% obsolete[,1],]
bp2 = bp2[order(bp2$p.adjust,decreasing=TRUE),]
#write.csv(bp2,"TRF_GO_CC_trimed_down_RNA_a.csv")

rk = character()
for(i in 1:(nrow(bp2)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp2)){
    pk = c(pk,sum(unlist(strsplit(bp2[i,8], "/")) %in% unlist(strsplit(bp2[j,8], "/")))/length(unlist(strsplit(bp2[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk,"ostavi")
bp2 = bp2[rk1 == "ostavi",]
#write.csv(bp2,"TRF_GO_CC_trimed_down_RNA_b.csv")

bp3 = ego3@result
#write.csv(bp3,"TRF_GO_MF_down_RNA.csv")

bp3 = bp3[bp3$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp3[,4]))
bp3 = bp3[(lk > 15) & (lk < 500),]
bp3 = bp3[bp3$Count > 4,]
bp3 = bp3[!bp3$ID %in% noliv[,1],]
bp3 = bp3[!bp3$ID %in% obsolete[,1],]
bp3 = bp3[order(bp3$p.adjust,decreasing=TRUE),]
#write.csv(bp3,"TRF_GO_MF_trimed_down_RNA_a.csv")

rk = character()
for(i in 1:(nrow(bp3)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp3)){
    pk = c(pk,sum(unlist(strsplit(bp3[i,8], "/")) %in% unlist(strsplit(bp3[j,8], "/")))/length(unlist(strsplit(bp3[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk,"ostavi")
bp3 = bp3[rk1 == "ostavi",]
#write.csv(bp3,"TRF_GO_MF_trimed_down_RNA_b.csv")

bp4 = kk1@result
for(i in 1:nrow(bp4)){
  bp4[i,8] <- paste(suppressMessages(bitr(unlist(strsplit(bp4[i,8], "/")),
                                          fromType="ENTREZID",
                                          toType="SYMBOL",
                                          OrgDb=org.Mm.eg.db)[,2]),
                    collapse="/"
)
}

#write.csv(bp4,"TRF_KEGG_down_RNA.csv")

bp4 = bp4[bp4$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp4[,4]))
bp4 = bp4[(lk > 15) & (lk < 500),]
bp4 = bp4[bp4$Count > 4,]
bp4 = bp4[!bp4$ID %in% noliv[,1],]
bp4 = bp4[!bp4$ID %in% obsolete[,1],]
bp4 = bp4[order(bp4$p.adjust,decreasing=TRUE),]

#write.csv(bp4,"TRF_KEGG_trimed_down_RNA_a.csv")

rk = character()
for(i in 1:(nrow(bp4)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp4)){
    pk = c(pk,sum(unlist(strsplit(bp4[i,8], "/")) %in% unlist(strsplit(bp4[j,8], "/")))/length(unlist(strsplit(bp4[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk,"ostavi")
bp4 = bp4[rk1 == "ostavi",]
#write.csv(bp4,"TRF_KEGG_trimed_down_RNA_b.csv")

###############################################################################
#
#   VSG study
#
###############################################################################

vsgdat = expression.matrix[,c(75,77:81,83:89)]

id = colnames(vsgdat)
condition = as.factor(c(rep("VSGHFD",6),rep("PFHFD",7)))
metaData = data.frame(id, condition)
dds <- DESeqDataSetFromMatrix(countData=vsgdat,
                              colData=metaData,
                              design=~condition
        )
featureData = data.frame(gene=rownames(dds))
mcols(dds) = DataFrame(mcols(dds), featureData)

design(dds) = ~ condition  
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds, 
                           fit='parametric'
)

keep = rowSums(counts(dds))  >= 10
dds = dds[keep,]

dds <- DESeq(dds)
res_vsghfd_vs_pfhfd = results(dds)

res = res_vsghfd_vs_pfhfd[order(res_vsghfd_vs_pfhfd$pvalue),]

resultsVSG = as.data.frame(dplyr::mutate(as.data.frame(res),
                                       significant=ifelse(res$pvalue < 0.05, "pvalue<0.05", "not signif.")), 
                         row.names=rownames(res))
resultsVSG = resultsVSG[(!is.na(resultsVSG$pvalue)),]

p <- ggplot2::ggplot(resultsVSG, ggplot2::aes(log2FoldChange, -log10(pvalue))) +
     ggplot2::geom_point(ggplot2::aes(col=significant)) + 
     ggplot2::scale_color_manual(values=c("grey", "forestgreen")) +
     guides(col=guide_legend(nrow=2)) + 
     ggplot2::ggtitle("") + 
     xlim(-5,5) + 
     ylim(0,20)
p

vsd <- assay(vst(dds,
                 fitType="local",
                 blind=TRUE)
)

df_pca = prcomp(t(vsd))
df_out = as.data.frame(df_pca$x)

p <- ggplot(df_out,aes(x=PC1,y=PC2,color=condition)) +
     ggplot2::scale_color_manual(values=c("green4", "green1"))
     geom_point(size=5)
p

#######################################
# UP
#######################################

sig.gene <- bitr(rownames(resultsVSG[which(resultsVSG$sig == "FDR<0.05" & resultsVSG$log2FoldChange > 0),]),
                 fromType="ENSEMBL",
                 toType="ENTREZID",
                 OrgDb=org.Mm.eg.db
)
kk1 <- enrichKEGG(gene=sig.gene[,2],
                  organism='mmu')
ego1 <- enrichGO(gene=sig.gene[,2],
                 OrgDb=org.Mm.eg.db,
                 ont="BP",
                 pAdjustMethod="BH",
                 readable=TRUE
)
ego2 <- enrichGO(gene=sig.gene[,2],
                 OrgDb=org.Mm.eg.db, 
                 ont="CC",
                 pAdjustMethod="BH", 
                 readable=TRUE
)
ego3 <- enrichGO(gene=sig.gene[,2],
                 OrgDb=org.Mm.eg.db,
                 ont="MF",
                 pAdjustMethod="BH", 
                 readable=TRUE
)
bp1 = ego1@result
setwd('/Users/viktorian.miok/Documents/consultation/Katarina/Integration_Data/RNA_seq_data/results/VSG')
#write.csv(bp1,"VSG_GO_BP_up_RNA.csv")

bp1 = bp1[bp1$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp1[,4]))
bp1 = bp1[(lk > 15) & (lk < 500),]
bp1 = bp1[bp1$Count > 4,]
bp1 = bp1[!bp1$ID %in% noliv[,1],]
bp1 = bp1[!bp1$ID %in% obsolete[,1],]
bp1 = bp1[order(bp1$p.adjust,decreasing=TRUE),]
#write.csv(bp1,"VSG_GO_BP_trimed_up_RNA_a.csv")

rk = character()
for(i in 1:(nrow(bp1)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp1)){
    pk = c(pk,sum(unlist(strsplit(bp1[i,8], "/")) %in% unlist(strsplit(bp1[j,8], "/")))/length(unlist(strsplit(bp1[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk,"ostavi")
bp1 = bp1[rk1 == "ostavi",]
#write.csv(bp1,"VSG_GO_BP_trimed_up_RNA_b.csv")

bp2 = ego2@result
#write.csv(bp2,"VSG_GO_CC_up_RNA.csv")

bp2 = bp2[bp2$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp2[,4]))
bp2 = bp2[(lk > 15) & (lk < 500),]
bp2 = bp2[bp2$Count > 4,]
bp2 = bp2[!bp2$ID %in% noliv[,1],]
bp2 = bp2[!bp2$ID %in% obsolete[,1],]
bp2 = bp2[order(bp2$p.adjust,decreasing=TRUE),]
#write.csv(bp2,"VSG_GO_CC_trimed_up_RNA_a.csv")

rk = character()
for(i in 1:(nrow(bp2)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp2)){
    pk = c(pk,sum(unlist(strsplit(bp2[i,8], "/")) %in% unlist(strsplit(bp2[j,8], "/")))/length(unlist(strsplit(bp2[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk,"ostavi")
bp2 = bp2[rk1 == "ostavi",]
#write.csv(bp2,"VSG_GO_CC_trimed_up_RNA_b.csv")

bp3 = ego3@result
#write.csv(bp3,"VSG_GO_MF_up_RNA.csv")

bp3 = bp3[bp3$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp3[,4]))
bp3 = bp3[(lk > 15) & (lk < 500),]
bp3 = bp3[bp3$Count > 4,]
bp3 = bp3[!bp3$ID %in% noliv[,1],]
bp3 = bp3[!bp3$ID %in% obsolete[,1],]
bp3 = bp3[order(bp3$p.adjust,decreasing=TRUE),]
#write.csv(bp3,"VSG_GO_MF_trimed_up_RNA_a.csv")

rk = character()
for(i in 1:(nrow(bp3)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp3)){
    pk = c(pk,sum(unlist(strsplit(bp3[i,8], "/")) %in% unlist(strsplit(bp3[j,8], "/")))/length(unlist(strsplit(bp3[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk,"ostavi")
bp3 = bp3[rk1 == "ostavi",]
#write.csv(bp3,"VSG_GO_MF_trimed_up_RNA_b.csv")

bp4 = kk1@result
for(i in 1:nrow(bp4)){
  bp4[i,8] <- paste(suppressMessages(bitr(unlist(strsplit(bp4[i,8], "/")),
                                          fromType="ENTREZID",
                                          toType="SYMBOL",
                                          OrgDb=org.Mm.eg.db)[,2]),
                    collapse="/"
)
}

#write.csv(bp4, "VSG_KEGG_up_RNA.csv")

bp4 = bp4[bp4$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp4[,4]))
bp4 = bp4[(lk > 15) & (lk < 500),]
bp4 = bp4[bp4$Count > 4,]
bp4 = bp4[!bp4$ID %in% noliv[,1],]
bp4 = bp4[!bp4$ID %in% obsolete[,1],]
bp4 = bp4[order(bp4$p.adjust,decreasing=TRUE),]

#write.csv(bp4, "VSG_KEGG_trimed_up_RNA_a.csv")

rk = character()
for(i in 1:(nrow(bp4)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp4)){
    pk = c(pk,sum(unlist(strsplit(bp4[i,8], "/")) %in% unlist(strsplit(bp4[j,8], "/")))/length(unlist(strsplit(bp4[i,8], "/"))))
  }
  if (sum(pk  >  0.4) > 0) {
    rk[i]<-"brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp4 = bp4[rk1 == "ostavi",]
#write.csv(bp4, "VSG_KEGG_trimed_up_RNA_b.csv")

#######################################
# DOWN
#######################################
sig.gene <- bitr(rownames(resultsVSG[which(resultsVSG$sig == "FDR<0.05" & resultsVSG$log2FoldChange < 0),]),
                 fromType="ENSEMBL",
                 toType="ENTREZID",
                 OrgDb=org.Mm.eg.db
)
kk1 <- enrichKEGG(gene=sig.gene[,2], 
                  organism='mmu'
)
ego1 <- enrichGO(gene=sig.gene[,2],
                 OrgDb=org.Mm.eg.db,
                 ont="BP",
                 pAdjustMethod="BH", 
                 readable=TRUE
)
ego2 <- enrichGO(gene=sig.gene[,2],
                 OrgDb=org.Mm.eg.db, 
                 ont="CC",
                 pAdjustMethod="BH",
                 readable=TRUE
)
ego3 <- enrichGO(gene=sig.gene[,2],
                 OrgDb=org.Mm.eg.db,
                 ont="MF", 
                 pAdjustMethod="BH", 
                 readable=TRUE
)
bp1 = ego1@result
setwd('/Users/viktorian.miok/Documents/consultation/Katarina/Integration_Data/RNA_seq_data/results/VSG')
#write.csv(bp1, "VSG_GO_BP_down_RNA.csv")

bp1 = bp1[bp1$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp1[,4]))
bp1 = bp1[(lk > 15) & (lk < 500),]
bp1 = bp1[bp1$Count > 4,]
bp1 = bp1[!bp1$ID %in% noliv[,1],]
bp1 = bp1[!bp1$ID %in% obsolete[,1],]
bp1 = bp1[order(bp1$p.adjust,decreasing=TRUE),]
#write.csv(bp1, "VSG_GO_BP_trimed_down_RNA_a.csv")

rk = character()
for(i in 1:(nrow(bp1)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp1)){
    pk = c(pk,sum(unlist(strsplit(bp1[i,8], "/")) %in% unlist(strsplit(bp1[j,8], "/")))/length(unlist(strsplit(bp1[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp1 = bp1[rk1 == "ostavi",]
#write.csv(bp1, "VSG_GO_BP_trimed_down_RNA_b.csv")

bp2 = ego2@result
#write.csv(bp2, "VSG_GO_CC_down_RNA.csv")

bp2 = bp2[bp2$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp2[,4]))
bp2 = bp2[(lk > 15) & (lk < 500),]
bp2 = bp2[bp2$Count > 4,]
bp2 = bp2[!bp2$ID %in% noliv[,1],]
bp2 = bp2[!bp2$ID %in% obsolete[,1],]
bp2 = bp2[order(bp2$p.adjust,decreasing=TRUE),]
#write.csv(bp2, "VSG_GO_CC_trimed_down_RNA_a.csv")

rk = character()
for(i in 1:(nrow(bp2)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp2)){
    pk = c(pk,sum(unlist(strsplit(bp2[i,8], "/")) %in% unlist(strsplit(bp2[j,8], "/")))/length(unlist(strsplit(bp2[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp2 = bp2[rk1 == "ostavi",]
#write.csv(bp2, "VSG_GO_CC_trimed_down_RNA_b.csv")

bp3 = ego3@result
#write.csv(bp3, "VSG_GO_MF_down_RNA.csv")

bp3 = bp3[bp3$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp3[,4]))
bp3 = bp3[(lk > 15) & (lk < 500),]
bp3 = bp3[bp3$Count > 4,]
bp3 = bp3[!bp3$ID %in% noliv[,1],]
bp3 = bp3[!bp3$ID %in% obsolete[,1],]
bp3 = bp3[order(bp3$p.adjust,decreasing=TRUE),]
#write.csv(bp3, "VSG_GO_MF_trimed_down_RNA_a.csv")

rk = character()
for(i in 1:(nrow(bp3)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp3)){
    pk = c(pk,sum(unlist(strsplit(bp3[i,8], "/")) %in% unlist(strsplit(bp3[j,8], "/")))/length(unlist(strsplit(bp3[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp3 = bp3[rk1 == "ostavi",]
#write.csv(bp3, "VSG_GO_MF_trimed_down_RNA_b.csv")

bp4 = kk1@result
for(i in 1:nrow(bp4)){
  bp4[i,8] <- paste(suppressMessages(bitr(unlist(strsplit(bp4[i,8], "/")),
                                          fromType="ENTREZID",
                                          toType="SYMBOL",
                                          OrgDb=org.Mm.eg.db)[,2]),
                    collapse="/"
)
}

#write.csv(bp4, "VSG_KEGG_down_RNA.csv")

bp4 = bp4[bp4$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp4[,4]))
bp4 = bp4[(lk > 15) & (lk < 500),]
bp4 = bp4[bp4$Count > 4,]
bp4 = bp4[!bp4$ID %in% noliv[,1],]
bp4 = bp4[!bp4$ID %in% obsolete[,1],]
bp4 = bp4[order(bp4$p.adjust,decreasing=TRUE),]

#write.csv(bp4, "VSG_KEGG_trimed_down_RNA_a.csv")

rk = character()
for(i in 1:(nrow(bp4)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp4)){
    pk = c(pk,sum(unlist(strsplit(bp4[i,8], "/")) %in% unlist(strsplit(bp4[j,8], "/")))/length(unlist(strsplit(bp4[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp4 = bp4[rk1 == "ostavi",]
#write.csv(bp4, "VSG_KEGG_trimed_down_RNA_b.csv")


###############################################################################
#
#   IF study
#
###############################################################################

ifdat = expression.matrix[,c(24,26:31,16:23)]  

id = colnames(ifdat)
condition = as.factor(c(rep("IFHFD",7),rep("ALHFD",8)))
metaData = data.frame(id, condition)
dds <- DESeqDataSetFromMatrix(countData=ifdat,
                              colData=metaData, 
                              design=~condition
)
featureData = data.frame(gene=rownames(dds))
mcols(dds) = DataFrame(mcols(dds), featureData)

design(dds) = ~condition  
dds = estimateSizeFactors(dds)
dds <- estimateDispersions(dds, 
                        fit='parametric'
)

keep = rowSums(counts(dds)) >= 10
dds = dds[keep,]

dds <- DESeq(dds)

res_ifhfd_vs_alhfd = results(dds, contrast=c("condition", "IFHFD", "ALHFD"))
res = res_ifhfd_vs_alhfd[order(res_ifhfd_vs_alhfd$pvalue),]
resultsIF = as.data.frame(dplyr::mutate(as.data.frame(res),
                                        significant=ifelse(res$pvalue < 0.05, "pvalue<0.05", "not signif.")),
                          row.names=rownames(res))
resultsIF = resultsIF[(!is.na(resultsIF$pvalue)),]

p <- ggplot2::ggplot(resultsIF, ggplot2::aes(log2FoldChange, -log10(pvalue))) +
     ggplot2::geom_point(ggplot2::aes(col=significant)) +
     ggplot2::scale_color_manual(values=c("grey", "firebrick")) +
     guides(col=guide_legend(nrow=2)) +
     ggplot2::ggtitle("") + 
     xlim(-5,5) + 
     ylim(0,20)
p

vsd <- assay(vst(dds,
                 fitType="local",
                 blind=TRUE)
)

df_pca <- prcomp(t(vsd))
df_out = as.data.frame(df_pca$x)

p <- ggplot(df_out,aes(x=PC1,y=PC2,color=condition)) +
     ggplot2::scale_color_manual(values=c("firebrick4", "firebrick1")) +
     geom_point(size=5)
p
#######################################
# UP
#######################################

sig.gene <- bitr(rownames(resultsIF[which(resultsIF$sig == "FDR<0.05" & resultsIF$log2FoldChange > 0),]),
                 fromType="ENSEMBL",
                 toType="ENTREZID",
                 OrgDb=org.Mm.eg.db
)
kk1 <- enrichKEGG(gene=sig.gene[,2],
                  organism='mmu'
)
ego1 <- enrichGO(gene=sig.gene[,2],
                 OrgDb=org.Mm.eg.db, 
                 ont="BP",
                 pAdjustMethod="BH", 
                 readable=TRUE
)
ego2 <- enrichGO(gene=sig.gene[,2],
                 OrgDb=org.Mm.eg.db, 
                 ont="CC", 
                 pAdjustMethod="BH",
                 readable=TRUE
)
ego3 <- enrichGO(gene=sig.gene[,2],
                 OrgDb=org.Mm.eg.db, 
                 ont="MF",
                 pAdjustMethod="BH",
                 readable=TRUE
)

bp1 = ego1@result
setwd('/Users/viktorian.miok/Documents/consultation/Katarina/Integration_Data/RNA_seq_data/results/IF')
#write.csv(bp1, "IF_GO_BP_up_RNA.csv")

bp1 = bp1[bp1$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp1[,4]))
bp1 = bp1[(lk > 15) & (lk < 500),]
bp1 = bp1[bp1$Count > 4,]
bp1 = bp1[!bp1$ID %in% noliv[,1],]
bp1 = bp1[!bp1$ID %in% obsolete[,1],]
bp1 = bp1[order(bp1$p.adjust,decreasing=TRUE),]
#write.csv(bp1, "IF_GO_BP_trimed_up_RNA_a.csv")

rk = character()
for(i in 1:(nrow(bp1)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp1)){
    pk = c(pk,sum(unlist(strsplit(bp1[i,8], "/")) %in% unlist(strsplit(bp1[j,8], "/")))/length(unlist(strsplit(bp1[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp1 = bp1[rk1 == "ostavi",]
#write.csv(bp1, "IF_GO_BP_trimed_up_RNA_b.csv")

bp2 = ego2@result
#write.csv(bp2, "IF_GO_CC_up_RNA.csv")

bp2 = bp2[bp2$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp2[,4]))
bp2 = bp2[(lk > 15) & (lk < 500),]
bp2 = bp2[bp2$Count > 4,]
bp2 = bp2[!bp2$ID %in% noliv[,1],]
bp2 = bp2[!bp2$ID %in% obsolete[,1],]
bp2 = bp2[order(bp2$p.adjust,decreasing=TRUE),]
#write.csv(bp2, "IF_GO_CC_trimed_up_RNA_a.csv")

rk = character()
for(i in 1:(nrow(bp2)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp2)){
    pk = c(pk,sum(unlist(strsplit(bp2[i,8], "/")) %in% unlist(strsplit(bp2[j,8], "/")))/length(unlist(strsplit(bp2[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp2 = bp2[rk1 == "ostavi",]
#write.csv(bp2, "IF_GO_CC_trimed_up_RNA_b.csv")

bp3 = ego3@result
#write.csv(bp3, "IF_GO_MF_up_RNA.csv")

bp3 = bp3[bp3$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp3[,4]))
bp3 = bp3[(lk > 15) & (lk < 500),]
bp3 = bp3[bp3$Count > 4,]
bp3 = bp3[!bp3$ID %in% noliv[,1],]
bp3 = bp3[!bp3$ID %in% obsolete[,1],]
bp3 = bp3[order(bp3$p.adjust,decreasing=TRUE),]
#write.csv(bp3, "IF_GO_MF_trimed_up_RNA_a.csv")

rk = character()
for(i in 1:(nrow(bp3)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp3)){
    pk = c(pk,sum(unlist(strsplit(bp3[i,8], "/")) %in% unlist(strsplit(bp3[j,8], "/")))/length(unlist(strsplit(bp3[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp3 = bp3[rk1 == "ostavi",]
#write.csv(bp3, "IF_GO_MF_trimed_up_RNA_b.csv")

bp4 = kk1@result
for(i in 1:nrow(bp4)){
  bp4[i,8] <- paste(suppressMessages(bitr(unlist(strsplit(bp4[i,8], "/")),
                                          fromType="ENTREZID",
                                          toType="SYMBOL",
                                          OrgDb=org.Mm.eg.db)[,2]),
                    collapse="/"
)
}

#write.csv(bp4, "IF_KEGG_up_RNA.csv")

bp4 = bp4[bp4$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp4[,4]))
bp4 = bp4[(lk > 15) & (lk < 500),]
bp4 = bp4[bp4$Count > 4,]
bp4 = bp4[!bp4$ID %in% noliv[,1],]
bp4 = bp4[!bp4$ID %in% obsolete[,1],]
bp4 = bp4[order(bp4$p.adjust,decreasing=TRUE),]

#write.csv(bp4, "IF_KEGG_trimed_up_RNA_a.csv")

rk = character()
for(i in 1:(nrow(bp4)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp4)){
    pk = c(pk,sum(unlist(strsplit(bp4[i,8], "/")) %in% unlist(strsplit(bp4[j,8], "/")))/length(unlist(strsplit(bp4[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp4 = bp4[rk1 == "ostavi",]
#write.csv(bp4, "IF_KEGG_trimed_up_RNA_b.csv")
#######################################
# DOWN
#######################################

sig.gene <- bitr(rownames(resultsIF[which(resultsIF$sig == "FDR<0.05" & resultsIF$log2FoldChange < 0),]),
                 fromType="ENSEMBL",
                 toType="ENTREZID",
                 OrgDb=org.Mm.eg.db
)
kk1 <- enrichKEGG(gene=sig.gene[,2],
                  organism='mmu'
)
ego1 <- enrichGO(gene=sig.gene[,2],
                 OrgDb=org.Mm.eg.db, 
                 ont="BP",
                 pAdjustMethod="BH",
                 readable=TRUE
)
ego2 <- enrichGO(gene=sig.gene[,2],
                 OrgDb=org.Mm.eg.db,
                 ont="CC", 
                 pAdjustMethod="BH", 
                 readable=TRUE
)
ego3 <- enrichGO(gene=sig.gene[,2], 
                 OrgDb=org.Mm.eg.db,
                 ont="MF",
                 pAdjustMethod="BH", 
                 readable=TRUE
)
bp1 = ego1@result
setwd('/Users/viktorian.miok/Documents/consultation/Katarina/Integration_Data/RNA_seq_data/results/IF')
#write.csv(bp1, "IF_GO_BP_down_RNA.csv")

bp1 = bp1[bp1$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp1[,4]))
bp1 = bp1[(lk > 15) & (lk < 500),]
bp1 = bp1[bp1$Count > 4,]
bp1 = bp1[!bp1$ID %in% noliv[,1],]
bp1 = bp1[!bp1$ID %in% obsolete[,1],]
bp1 = bp1[order(bp1$p.adjust,decreasing=TRUE),]# (as.numeric(sub("\\/.*", "", bp1[,4]))
#write.csv(bp1, "IF_GO_BP_trimed_down_RNA_a.csv")

rk = character()
for(i in 1:(nrow(bp1)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp1)){
    pk = c(pk,sum(unlist(strsplit(bp1[i,8], "/")) %in% unlist(strsplit(bp1[j,8], "/")))/length(unlist(strsplit(bp1[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp1 = bp1[rk1 == "ostavi",]
#write.csv(bp1, "IF_GO_BP_trimed_down_RNA_b.csv")

bp2 = ego2@result
#write.csv(bp2, "IF_GO_CC_down_RNA.csv")

bp2 = bp2[bp2$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp2[,4]))
bp2 = bp2[(lk > 15) & (lk < 500),]
bp2 = bp2[bp2$Count > 4,]
bp2 = bp2[!bp2$ID %in% noliv[,1],]
bp2 = bp2[!bp2$ID %in% obsolete[,1],]
bp2 = bp2[order(bp2$p.adjust,decreasing=TRUE),]
#write.csv(bp2, "IF_GO_CC_trimed_down_RNA_a.csv")

rk = character()
for(i in 1:(nrow(bp2)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp2)){
    pk = c(pk,sum(unlist(strsplit(bp2[i,8], "/")) %in% unlist(strsplit(bp2[j,8], "/")))/length(unlist(strsplit(bp2[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp2 = bp2[rk1 == "ostavi",]
#write.csv(bp2, "IF_GO_CC_trimed_down_RNA_b.csv")

bp3 = ego3@result
#write.csv(bp3, "IF_GO_MF_down_RNA.csv")

bp3 = bp3[bp3$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp3[,4]))
bp3 = bp3[(lk > 15) & (lk < 500),]
bp3 = bp3[bp3$Count > 4,]
bp3 = bp3[!bp3$ID %in% noliv[,1],]
bp3 = bp3[!bp3$ID %in% obsolete[,1],]
bp3 = bp3[order(bp3$p.adjust,decreasing=TRUE),]
#write.csv(bp3, "IF_GO_MF_trimed_down_RNA_a.csv")

rk = character()
for(i in 1:(nrow(bp3)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp3)){
    pk = c(pk,sum(unlist(strsplit(bp3[i,8], "/")) %in% unlist(strsplit(bp3[j,8], "/")))/length(unlist(strsplit(bp3[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp3 = bp3[rk1 == "ostavi",]
#write.csv(bp3, "IF_GO_MF_trimed_down_RNA_b.csv")

bp4 = kk1@result
for(i in 1:nrow(bp4)){
  bp4[i,8] <- paste(suppressMessages(bitr(unlist(strsplit(bp4[i,8], "/")),
                                          fromType="ENTREZID",
                                          toType="SYMBOL",
                                          OrgDb=org.Mm.eg.db)[,2]),
                    collapse="/"
)
}

#write.csv(bp4, "IF_KEGG_down_RNA.csv")

bp4 = bp4[bp4$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp4[,4]))
bp4 = bp4[(lk > 15) & (lk < 500),]
bp4 = bp4[bp4$Count > 4,]
bp4 = bp4[!bp4$ID %in% noliv[,1],]
bp4 = bp4[!bp4$ID %in% obsolete[,1],]
bp4 = bp4[order(bp4$p.adjust,decreasing=TRUE),]

#write.csv(bp4, "IF_KEGG_trimed_down_RNA_a.csv")

rk = character()
for(i in 1:(nrow(bp4)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp4)){
    pk = c(pk,sum(unlist(strsplit(bp4[i,8], "/")) %in% unlist(strsplit(bp4[j,8], "/")))/length(unlist(strsplit(bp4[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp4 = bp4[rk1 == "ostavi",]
#write.csv(bp4, "IF_KEGG_trimed_down_RNA_b.csv")

###############################################################################
#
#     Integrated analysis
#
###############################################################################

setwd('/Users/viktorian.miok/Documents/consultation/Katarina/Integration_Data/RNA_seq_data')
load("ProteomOutput.RData")

protVSG = ProteomOutput$diffVSG
protTRF = ProteomOutput$diffTRF
protIF = ProteomOutput$diffIF

#######################################
# IF study
#######################################
IF = cbind(rownames(resultsIF), resultsIF[,c(2,5,6)])
colnames(IF) = c("ENSEMBL", "log2FoldChange", "pvalue", "padj")
IF[,1] = as.character(IF[,1])

sig.gene <- bitr(IF[,1],
                 fromType="ENSEMBL",
                 toType="SYMBOL",
                 OrgDb=org.Mm.eg.db
)
colnames(sig.gene) = c("ENSEMBL", "GeneName")

rnaIF <- merge(IF,
               sig.gene,
               by="ENSEMBL")[,c(5,2,3)] # ,4

#######################################
# VSG study
#######################################
VSG = cbind(rownames(resultsVSG), resultsVSG[,c(2,5,6)])
colnames(VSG) = c("ENSEMBL", "log2FoldChange", "pvalue", "padj")
VSG[,1] = as.character(VSG[,1])

sig.gene <- bitr(VSG[,1],
                 fromType="ENSEMBL",
                 toType="SYMBOL",
                 OrgDb=org.Mm.eg.db
)
colnames(sig.gene)=c("ENSEMBL", "GeneName")

rnaVSG <- merge(VSG,
                sig.gene,
                by="ENSEMBL")[,c(5,2,3)] # ,4

#######################################
# TRF study
#######################################

TRF = cbind(rownames(resultsTRF), resultsTRF[,c(2,5,6)])
colnames(TRF) = c("ENSEMBL", "log2FoldChange", "pvalue", "padj")
TRF[,1] = as.character(TRF[,1])

sig.gene <- bitr(TRF[,1],
                 fromType="ENSEMBL",
                 toType="SYMBOL",
                 OrgDb=org.Mm.eg.db
)
colnames(sig.gene) = c("ENSEMBL", "GeneName")

rnaTRF <- merge(TRF,
                sig.gene,
                by="ENSEMBL")[,c(5,2,3)] #,4


VSGfin <- merge(protVSG, 
                rnaVSG,
                by="GeneName"
)
TRFfin <- merge(protTRF,
                rnaTRF,
                by="GeneName"
)
IFfin <- merge(protIF,
               rnaIF,
               by="GeneName"
)

VSGfin["log10Pval_Trans"] = -log10(VSGfin$pvalue.y)*sign(VSGfin$log2FoldChange.y)
VSGfin["log10Pval_Prote"] = -log10(VSGfin$pvalue.x)*sign(VSGfin$log2FoldChange.x)

TRFfin["log10Pval_Trans"] = -log10(TRFfin$pvalue.y)*sign(TRFfin$log2FoldChange.y)
TRFfin["log10Pval_Prote"] = -log10(TRFfin$pvalue.x)*sign(TRFfin$log2FoldChange.x)

IFfin["log10Pval_Trans"] = -log10(IFfin$pvalue.y)*sign(IFfin$log2FoldChange.y)
IFfin["log10Pval_Prote"] = -log10(IFfin$pvalue.x)*sign(IFfin$log2FoldChange.x)

cor(VSGfin$log10Pval_Trans,
    VSGfin$log10Pval_Prote,
    method="spearman"
)
cor(TRFfin$log10Pval_Trans, 
    TRFfin$log10Pval_Prote, 
    method="spearman")
cor(IFfin$log10Pval_Trans, 
    IFfin$log10Pval_Prote, 
    method="spearman")

cutoff = 0.05

protVSG1 = protVSG[which(protVSG$pvalue < cutoff),]
protTRF1 = protTRF[which(protTRF$pvalue < cutoff),]
protIF1 = protIF[which(protIF$pvalue < cutoff),]

rnaVSG1 = rnaVSG[which(rnaVSG$pvalue < cutoff),]
rnaTRF1 = rnaTRF[which(rnaTRF$pvalue < cutoff),]
rnaIF1 = rnaIF[which(rnaIF$pvalue < cutoff),]

VSGfin1 <- merge(protVSG1,
                 rnaVSG1, 
                 by="GeneName"
)
TRFfin1 <- merge(protTRF1,
                 rnaTRF1, 
                 by="GeneName"
)
IFfin1 <- merge(protIF1, 
                rnaIF1,
                by="GeneName"
)
VSGfin1["log10Pval_Trans"] = -log10(VSGfin1$pvalue.y)*sign(VSGfin1$log2FoldChange.y)
VSGfin1["log10Pval_Prote"] = -log10(VSGfin1$pvalue.x)*sign(VSGfin1$log2FoldChange.x)

TRFfin1["log10Pval_Trans"] = -log10(TRFfin1$pvalue.y)*sign(TRFfin1$log2FoldChange.y)
TRFfin1["log10Pval_Prote"] = -log10(TRFfin1$pvalue.x)*sign(TRFfin1$log2FoldChange.x)

IFfin1["log10Pval_Trans"] = -log10(IFfin1$pvalue.y)*sign(IFfin1$log2FoldChange.y)
IFfin1["log10Pval_Prote"] = -log10(IFfin1$pvalue.x)*sign(IFfin1$log2FoldChange.x)
#############################################################################################################################
#
#       Integrate IF and TRF
#
############################################################################################################################
TRF_IFrna <- merge(rnaTRF,
                   rnaIF,
                   by="GeneName"
)
TRF_IFprot <- merge(protTRF, 
                    protIF, 
                    by="GeneName"
)
TRF_IFrna["log10Pval_TRF"] = -log10(TRF_IFrna$pvalue.y)*sign(TRF_IFrna$log2FoldChange.y)
TRF_IFrna["log10Pval_IF"] = -log10(TRF_IFrna$pvalue.x)*sign(TRF_IFrna$log2FoldChange.x)

TRF_IFprot["log10Pval_TRF"] = -log10(TRF_IFprot$pvalue.y)*sign(TRF_IFprot$log2FoldChange.y)
TRF_IFprot["log10Pval_IF"] = -log10(TRF_IFprot$pvalue.x)*sign(TRF_IFprot$log2FoldChange.x)

cor(TRF_IFrna$log10Pval_TRF, 
    TRF_IFrna$log10Pval_IF,
    method="spearman"
) # -0.25
cor(TRF_IFprot$log10Pval_TRF,
    TRF_IFprot$log10Pval_IF,
    method="spearman"
) # 0.44
###############################################################################
#       plot transcriptomics IF - TRF 
###############################################################################

TRF_IFrna_gen = c("Nfil3","Nr1d2","Nr1d1","Per3","Per2","Rorc","Cry1","Arg1","Asl","Hgd","Fah","Aldh4a1","Prodh","Abcg8","Npc1")
TRF_IFrna$significance <- ifelse((TRF_IFrna$log10Pval_IF > 1.30103) & (TRF_IFrna$log10Pval_TRF > 1.30103), "commonly", 
                                 ifelse((TRF_IFrna$log10Pval_IF > 1.30103) & (TRF_IFrna$log10Pval_TRF < -1.30103), "opposite", 
                                        ifelse((TRF_IFrna$log10Pval_IF < -1.30103) & (TRF_IFrna$log10Pval_TRF > 1.30103), "opposite",
                                              ifelse((TRF_IFrna$log10Pval_IF < -1.30103) & (TRF_IFrna$log10Pval_TRF < -1.30103), "commonly", "non-significant")
                                        )      
                                 )
)
  
#  ifelse((TRF_IFrna$pvalue.x < 0.05) & (TRF_IFrna$pvalue.y < 0.05), "p-value<0.05", "p-value>0.05")
TRF_IFrna$label <- ifelse(TRF_IFrna$GeneName %in% TRF_IFrna_gen,
                          "name",
                          ""
)
for(i in 1:nrow(TRF_IFrna)){
  if(TRF_IFrna$label[i] == "name") TRF_IFrna$label[i] = TRF_IFrna$GeneName[i]
}

setwd("/Users/viktorian.miok/Documents/consultation/Katarina/Integration_Data/RNA_seq_data/results/TRF_IF/")
#pdf("TRF_IF_rna.pdf") 
b <- ggplot(TRF_IFrna, aes(x=log10Pval_IF, y=log10Pval_TRF)) +
           geom_point(aes(color=significance), show.legend=FALSE) +
           scale_color_manual(values=c("mediumpurple3", "grey", "goldenrod3")) + 
           geom_hline(yintercept=c(-log10(cutoff), log10(cutoff)), linetype="dashed") +
           geom_vline( xintercept=c(-log10(cutoff), log10(cutoff)), linetype="dashed") +
           xlab("-Log10Pval(IF)*sign(FC)") + 
           ylab("-Log10Pval(TRF)*sign(FC)")+
           annotate(geom="text",
                    x=9,
                    y=16,
                    label=expression(correlaton:rho  ==  -0.25), color="blue") + 
           ggtitle("IF-TRF transcriptomics study")+
           geom_text_repel(aes(label=label),
                        size=5,
                        box.padding=unit(1.3, "lines"),
                        point.padding=unit(0.1, "lines"),
                        segment.size=0.15) 
#dev.off()

########################################################
#   tables of significant transcriptomics IF - TRF  
########################################################
colnames(TRF_IFrna) = c("GeneName","log2FoldChange.TRF","pvalue.TRF","log2FoldChange.IF","pvalue.IF",
                        "-Log10Pval(TRF)*sign(FC)","-Log10Pval(IF)*sign(FC)","significance","label")
TRF_IFrnaIntegrated_common = TRF_IFrna[TRF_IFrna$significance == "commonly",1:7]
TRF_IFrnaIntegrated_opposite = TRF_IFrna[TRF_IFrna$significance == "opposite",1:7]
TRF_IFrnaIntegrated = TRF_IFrna[,1:7]

#setwd("/Users/viktorian.miok/Documents/consultation/Katarina/Integration_Data/RNA_seq_data/results/TRF_IF")
#write.csv(TRF_IFrnaIntegrated_common, "TRF_IFrnaIntegrated_common.csv")
#write.csv(TRF_IFrnaIntegrated_opposite, "TRF_IFrnaIntegrated_opposite.csv")
#write.csv(TRF_IFrnaIntegrated, "TRF_IFrnaIntegrated.csv")

###############################################################################
#       proteomics IF - TRF 
###############################################################################

TRF_IFprot_gen = c("Aars","Iars","Nars","Ndufa13","Ndufa9","Ndufb11","Ndufb3","Ndufb4","Apoa2","Apoe")
TRF_IFprot$significance=ifelse((TRF_IFprot$log10Pval_IF > 1.30103) & (TRF_IFprot$log10Pval_TRF > 1.30103), "commonly", 
                                 ifelse((TRF_IFprot$log10Pval_IF > 1.30103) & (TRF_IFprot$log10Pval_TRF < -1.30103), "opposite", 
                                        ifelse((TRF_IFprot$log10Pval_IF < -1.30103) & (TRF_IFprot$log10Pval_TRF > 1.30103), "opposite",
                                               ifelse((TRF_IFprot$log10Pval_IF < -1.30103) & (TRF_IFprot$log10Pval_TRF < -1.30103), "commonly", "non-significant")
                                        )      
                                 )
)
TRF_IFprot$label <- ifelse(TRF_IFprot$GeneName %in% TRF_IFprot_gen,
                           "name",
                           ""
)
for(i in 1:nrow(TRF_IFprot)){
  if(TRF_IFprot$label[i] == "name") TRF_IFprot$label[i] = TRF_IFprot$GeneName[i]
}

setwd("/Users/viktorian.miok/Documents/consultation/Katarina/Integration_Data/Proteomics_data/results/TRF_IF/")
#pdf("TRF_IF_protein.pdf") 
b <- ggplot(TRF_IFprot, aes(x=log10Pval_IF, y=log10Pval_TRF)) +
           geom_point(aes(color=significance), show.legend=FALSE) +
           scale_color_manual(values=c("mediumpurple3", "grey", "goldenrod3")) +
           geom_hline(yintercept=c(-log10(cutoff),log10(cutoff)), linetype="dashed") +
           geom_vline( xintercept=c(-log10(cutoff),log10(cutoff)), linetype="dashed") + 
           xlab("-Log10Pval(IF)*sign(FC)") + 
           ylab("-Log10Pval(TRF)*sign(FC)")+
           annotate(geom="text", x=-3.5, y=6, label=expression(correlaton:rho == 0.44), color="blue") +
           ggtitle("IF-TRF protein study") +
           geom_text_repel(aes(label=label),
                           size=5,
                           box.padding=unit(1.3, "lines"), 
                           point.padding=unit(0.1, "lines"),
                           segment.size=0.15) 
#dev.off()

########################################################
#   tables of significant proteomics IF - TRF  
########################################################
colnames(TRF_IFprot) = c("GeneName","log2FoldChange.TRF","pvalue.TRF","log2FoldChange.IF","pvalue.IF",
                         "-Log10Pval(TRF)*sign(FC)","-Log10Pval(IF)*sign(FC)","significance","label")
TRF_IFprotIntegrated_common = TRF_IFprot[TRF_IFprot$significance == "commonly", 1:7]
TRF_IFprotIntegrated_opposite = TRF_IFprot[TRF_IFprot$significance == "opposite", 1:7]
TRF_IFprotIntegrated = TRF_IFprot[,1:7]

#setwd("/Users/viktorian.miok/Documents/consultation/Katarina/Integration_Data/Proteomics_data/results/TRF_IF")
#write.csv(TRF_IFprotIntegrated_common, "TRF_IFproteinIntegrated_common.csv")
#write.csv(TRF_IFprotIntegrated_opposite, "TRF_IFproteinIntegrated_opposite.csv")
#write.csv(TRF_IFprotIntegrated, "TRF_IFproteinIntegrated.csv")

#############################################################################################################################
#
#       Integrate VSG and TRF
#
############################################################################################################################
VSG_TRFrna <- merge(rnaVSG, 
                    rnaTRF,
                    by="GeneName"
)
VSG_TRFprot <- merge(protVSG,
                     protTRF, 
                     by="GeneName"
)

VSG_TRFrna["log10Pval_VSG"] = -log10(VSG_TRFrna$pvalue.x)*sign(VSG_TRFrna$log2FoldChange.x)
VSG_TRFrna["log10Pval_TRF"] = -log10(VSG_TRFrna$pvalue.y)*sign(VSG_TRFrna$log2FoldChange.y)

VSG_TRFprot["log10Pval_VSG"] = -log10(VSG_TRFprot$pvalue.x)*sign(VSG_TRFprot$log2FoldChange.x)
VSG_TRFprot["log10Pval_TRF"] = -log10(VSG_TRFprot$pvalue.y)*sign(VSG_TRFprot$log2FoldChange.y)

cor(VSG_TRFrna$log10Pval_VSG,
    VSG_TRFrna$log10Pval_TRF,
    method="p"
) # -0.25
cor(VSG_TRFprot$log10Pval_VSG, 
    VSG_TRFprot$log10Pval_TRF, 
    method="p"
) # 0.44
###############################################################################
#       plot transcriptomics VSG - TRF 
###############################################################################

VSG_TRFrna_gen = c("Nfil3","Nr1d2","Nr1d1","Per3","Per2","Rorc","Cry1","Arg1","Asl","Hgd","Fah","Aldh4a1","Prodh","Abcg8","Npc1")
VSG_TRFrna$significance <- ifelse((VSG_TRFrna$log10Pval_VSG > 1.30103) & (VSG_TRFrna$log10Pval_TRF > 1.30103), "commonly", 
                                 ifelse((VSG_TRFrna$log10Pval_VSG > 1.30103) & (VSG_TRFrna$log10Pval_TRF < -1.30103), "opposite", 
                                        ifelse((VSG_TRFrna$log10Pval_VSG < -1.30103) & (VSG_TRFrna$log10Pval_TRF > 1.30103), "opposite",
                                               ifelse((VSG_TRFrna$log10Pval_VSG < -1.30103) & (VSG_TRFrna$log10Pval_TRF < -1.30103), "commonly", "non-significant")
                                        )      
                                 )
)

VSG_TRFrna$label = ifelse(VSG_TRFrna$GeneName %in% VSG_TRFrna_gen, "name", "")
for(i in 1:nrow(VSG_TRFrna)){
  if(VSG_TRFrna$label[i] == "name") VSG_TRFrna$label[i] = VSG_TRFrna$GeneName[i]
}

setwd("/Users/viktorian.miok/Documents/consultation/Katarina/Integration_Data/RNA_seq_data/results/VSG_TRF/")
pdf("VSG_TRF_rna.pdf") 
b <- ggplot(VSG_TRFrna, aes(x=log10Pval_VSG, y=log10Pval_TRF)) +
           geom_point(aes(color=significance), show.legend=FALSE) +
           scale_color_manual(values=c("mediumpurple3", "grey", "goldenrod3")) + 
           geom_hline(yintercept=c(-log10(cutoff),log10(cutoff)),linetype="dashed") +
           geom_vline( xintercept=c(-log10(cutoff),log10(cutoff)),linetype="dashed") +
           xlab("-Log10Pval(VSG)*sign(FC)") + 
           ylab("-Log10Pval(TRF)*sign(FC)") +
           annotate(geom="text", x=-5.5, y=16, label=expression(correlaton:rho == 0.17), color="blue") + 
           ggtitle("VSG-TRF transcriptomics study") +
           geom_text_repel(aes(label=label),
                           size=5,
                           box.padding=unit(0.5, "lines"),
                           point.padding=unit(0.1, "lines"),
                           segment.size=0.15) 
dev.off()

########################################################
#   tables of significant transcriptomics VSG - TRF  
########################################################
colnames(VSG_TRFrna) = c("GeneName","log2FoldChange.VSG","pvalue.VSG","log2FoldChange.TRF","pvalue.TRF",
                         "-Log10Pval(VSG)*sign(FC)","-Log10Pval(TRF)*sign(FC)","significance","label")
VSG_TRFrnaIntegrated_common = VSG_TRFrna[VSG_TRFrna$significance == "commonly", 1:7]
VSG_TRFrnaIntegrated_opposite = VSG_TRFrna[VSG_TRFrna$significance == "opposite", 1:7]
VSG_TRFrnaIntegrated = VSG_TRFrna[,1:7]

write.csv(VSG_TRFrnaIntegrated_common, "VSG_TRFrnaIntegrated_common.csv")
write.csv(VSG_TRFrnaIntegrated_opposite, "VSG_TRFrnaIntegrated_opposite.csv")
write.csv(VSG_TRFrnaIntegrated, "VSG_TRFrnaIntegrated.csv")

###############################################################################
#       proteomics VSG - TRF 
###############################################################################

VSG_TRFprot_gen=c("Aars","Iars","Nars","Ndufa13","Ndufa9","Ndufb11","Ndufb3","Ndufb4","Apoa2","Apoe")
VSG_TRFprot$significance <- ifelse((VSG_TRFprot$log10Pval_VSG > 1.30103) & (VSG_TRFprot$log10Pval_TRF > 1.30103), "commonly", 
                                  ifelse((VSG_TRFprot$log10Pval_VSG > 1.30103) & (VSG_TRFprot$log10Pval_TRF < -1.30103), "opposite", 
                                         ifelse((VSG_TRFprot$log10Pval_VSG < -1.30103) & (VSG_TRFprot$log10Pval_TRF > 1.30103), "opposite",
                                                ifelse((VSG_TRFprot$log10Pval_VSG < -1.30103) & (VSG_TRFprot$log10Pval_TRF < -1.30103), "commonly", "non-significant")
                                         )      
                                  )
)
VSG_TRFprot$label <- ifelse(VSG_TRFprot$GeneName %in% VSG_TRFprot_gen,
                            "name",
                            ""
)
for(i in 1:nrow(VSG_TRFprot)){
  if(VSG_TRFprot$label[i] == "name") VSG_TRFprot$label[i] = VSG_TRFprot$GeneName[i]
}

setwd("/Users/viktorian.miok/Documents/consultation/Katarina/Integration_Data/Proteomics_data/results/VSG_TRF/")
pdf("VSG_TRF_protein.pdf") 
b <- ggplot(VSG_TRFprot, aes(x=log10Pval_VSG, y=log10Pval_TRF)) + 
           geom_point(aes(color=significance), show.legend=FALSE) +
           scale_color_manual(values=c("mediumpurple3", "grey", "goldenrod3")) + 
           geom_hline(yintercept=c(-log10(cutoff), log10(cutoff)), linetype="dashed") +
           geom_vline( xintercept=c(-log10(cutoff), log10(cutoff)), linetype="dashed") +
           xlab("-Log10Pval(VSG)*sign(FC)") + 
           ylab("-Log10Pval(TRF)*sign(FC)")+
           annotate(geom="text", x=-4, y=6, label=expression(correlaton:rho == 0.3), color="blue") +
           ggtitle("VSG-TRF protein study")+
           geom_text_repel(aes(label=label),
                           size=5, 
                           box.padding=unit(1.3, "lines"), 
                           point.padding=unit(0.1, "lines"),
                           segment.size=0.15) 
dev.off()

########################################################
#   tables of significant proteomics VSG - TRF  
########################################################
colnames(VSG_TRFprot) = c("GeneName","log2FoldChange.VSG","pvalue.VSG","log2FoldChange.TRF","pvalue.TRF",
                          "-Log10Pval(VSG)*sign(FC)","-Log10Pval(TRF)*sign(FC)","significance","label")
VSG_TRFprotIntegrated_common = VSG_TRFprot[VSG_TRFprot$significance == "commonly", 1:7]
VSG_TRFprotIntegrated_opposite = VSG_TRFprot[VSG_TRFprot$significance == "opposite", 1:7]
VSG_TRFprotIntegrated=VSG_TRFprot[,1:7]

write.csv(VSG_TRFprotIntegrated_common, "VSG_TRFproteinIntegrated_common.csv")
write.csv(VSG_TRFprotIntegrated_opposite, "VSG_TRFproteinIntegrated_opposite.csv")
write.csv(VSG_TRFprotIntegrated, "VSG_TRFproteinIntegrated.csv")


#############################################################################################################################
#
#       Integrate VSG and IF
#
############################################################################################################################
VSG_IFrna <- merge(rnaVSG, 
                   rnaIF,
                   by="GeneName"
)
VSG_IFprot <- merge(protVSG,
                    protIF, 
                    by="GeneName"
)

VSG_IFrna["log10Pval_VSG"] = -log10(VSG_IFrna$pvalue.x)*sign(VSG_IFrna$log2FoldChange.x)
VSG_IFrna["log10Pval_IF"] = -log10(VSG_IFrna$pvalue.y)*sign(VSG_IFrna$log2FoldChange.y)

VSG_IFprot["log10Pval_VSG"] = -log10(VSG_IFprot$pvalue.x)*sign(VSG_IFprot$log2FoldChange.x)
VSG_IFprot["log10Pval_IF"] = -log10(VSG_IFprot$pvalue.y)*sign(VSG_IFprot$log2FoldChange.y)

cor(VSG_IFrna$log10Pval_VSG, 
    VSG_IFrna$log10Pval_IF, 
    method="p"
) # -0.25
cor(VSG_IFprot$log10Pval_VSG,
    VSG_IFprot$log10Pval_IF,
    method="p"
) # 0.44
###############################################################################
#       plot transcriptomics VSG - IF 
###############################################################################

VSG_IFrna_gen = c("Nfil3","Nr1d2","Nr1d1","Per3","Per2","Rorc","Cry1","Arg1","Asl","Hgd","Fah","Aldh4a1","Prodh","Abcg8","Npc1")
VSG_IFrna$significance <- ifelse((VSG_IFrna$log10Pval_VSG > 1.30103) & (VSG_IFrna$log10Pval_IF > 1.30103), "commonly", 
                                  ifelse((VSG_IFrna$log10Pval_VSG > 1.30103) & (VSG_IFrna$log10Pval_IF < -1.30103), "opposite", 
                                         ifelse((VSG_IFrna$log10Pval_VSG < -1.30103) & (VSG_IFrna$log10Pval_IF > 1.30103), "opposite",
                                                ifelse((VSG_IFrna$log10Pval_VSG < -1.30103) & (VSG_IFrna$log10Pval_IF < -1.30103), "commonly", "non-significant")
                                         )      
                                  )
)

VSG_IFrna$label <- ifelse(VSG_IFrna$GeneName %in% VSG_IFrna_gen,
                          "name",
                          ""
)
for(i in 1:nrow(VSG_IFrna)){
  if(VSG_IFrna$label[i] == "name") VSG_IFrna$label[i] = VSG_IFrna$GeneName[i]
}

setwd("/Users/viktorian.miok/Documents/consultation/Katarina/Integration_Data/RNA_seq_data/results/VSG_IF/")
pdf("VSG_IF_rna.pdf") 
b <- ggplot(VSG_IFrna, aes(x=log10Pval_VSG, y=log10Pval_IF)) + 
           geom_point(aes(color=significance), show.legend=FALSE) +
           scale_color_manual(values=c("mediumpurple3", "grey", "goldenrod3")) +
           geom_hline(yintercept=c(-log10(cutoff), log10(cutoff)), linetype="dashed") +
           geom_vline( xintercept=c(-log10(cutoff), log10(cutoff)), linetype="dashed") +
           xlab("-Log10Pval(VSG)*sign(FC)") +
           ylab("-Log10Pval(IF)*sign(FC)")+
           annotate(geom="text", x=-5.5, y=16, label=expression(correlaton:rho == -0.08), color="blue") + 
           ggtitle("VSG-IF transcriptomics study")+
           geom_text_repel(aes(label=label),
                           size=5,
                           box.padding=unit(0.5, "lines"),
                           point.padding=unit(0.1, "lines"),
                           segment.size=0.15) 
dev.off()

########################################################
#   tables of significant transcriptomics VSG - IF  
########################################################
colnames(VSG_IFrna)=c("GeneName","log2FoldChange.VSG","pvalue.VSG","log2FoldChange.IF","pvalue.IF",
                      "-Log10Pval(VSG)*sign(FC)","-Log10Pval(IF)*sign(FC)","significance","label")
VSG_IFrnaIntegrated_common = VSG_IFrna[VSG_IFrna$significance == "commonly", 1:7]
VSG_IFrnaIntegrated_opposite = VSG_IFrna[VSG_IFrna$significance == "opposite", 1:7]
VSG_IFrnaIntegrated = VSG_IFrna[,1:7]

write.csv(VSG_IFrnaIntegrated_common, "VSG_IFrnaIntegrated_common.csv")
write.csv(VSG_IFrnaIntegrated_opposite, "VSG_IFrnaIntegrated_opposite.csv")
write.csv(VSG_IFrnaIntegrated, "VSG_IFrnaIntegrated.csv")

###############################################################################
#       proteomics VSG - IF 
###############################################################################

VSG_IFprot_gen = c("Aars","Iars","Nars","Ndufa13","Ndufa9","Ndufb11","Ndufb3","Ndufb4","Apoa2","Apoe")
VSG_IFprot$significance <- ifelse((VSG_IFprot$log10Pval_VSG > 1.30103) & (VSG_IFprot$log10Pval_IF > 1.30103), "commonly", 
                                   ifelse((VSG_IFprot$log10Pval_VSG > 1.30103) & (VSG_IFprot$log10Pval_IF < -1.30103), "opposite", 
                                          ifelse((VSG_IFprot$log10Pval_VSG < -1.30103) & (VSG_IFprot$log10Pval_IF > 1.30103), "opposite",
                                                 ifelse((VSG_IFprot$log10Pval_VSG < -1.30103) & (VSG_IFprot$log10Pval_IF < -1.30103), "commonly", "non-significant")
                                          )      
                                   )
)
VSG_IFprot$label <- ifelse(VSG_IFprot$GeneName %in% VSG_IFprot_gen,
                           "name", 
                           ""
)
for(i in 1:nrow(VSG_IFprot)){
  if(VSG_IFprot$label[i] == "name") VSG_IFprot$label[i] = VSG_IFprot$GeneName[i]
}

setwd("/Users/viktorian.miok/Documents/consultation/Katarina/Integration_Data/Proteomics_data/results/VSG_IF/")
pdf("VSG_IF_protein.pdf") 
b <- ggplot(VSG_IFprot, aes(x=log10Pval_VSG, y=log10Pval_IF)) + 
           geom_point(aes(color=significance), show.legend=FALSE) +
           scale_color_manual(values=c("mediumpurple3", "grey", "goldenrod3")) +
           geom_hline(yintercept=c(-log10(cutoff),log10(cutoff)),linetype="dashed") +
           geom_vline( xintercept=c(-log10(cutoff),log10(cutoff)),linetype="dashed") +
           xlab("-Log10Pval(VSG)*sign(FC)") + 
           ylab("-Log10Pval(IF)*sign(FC)")+
           annotate(geom="text", x=-4, y=6, label=expression(correlaton:rho == -0.11), color="blue") + 
           ggtitle("VSG-IF protein study")+
           geom_text_repel(aes(label=label),
                           size=5,
                           box.padding=unit(1.3, "lines"),
                           point.padding=unit(0.1, "lines"),
                           segment.size=0.15) 
dev.off()

########################################################
#   tables of significant proteomics VSG - IF  
########################################################
colnames(VSG_IFprot) = c("GeneName","log2FoldChange.VSG","pvalue.VSG","log2FoldChange.IF","pvalue.IF",
                         "-Log10Pval(VSG)*sign(FC)","-Log10Pval(IF)*sign(FC)","significance","label")
VSG_IFprotIntegrated_common = VSG_IFprot[VSG_IFprot$significance == "commonly", 1:7]
VSG_IFprotIntegrated_opposite = VSG_IFprot[VSG_IFprot$significance == "opposite", 1:7]
VSG_IFprotIntegrated = VSG_IFprot[,1:7]

write.csv(VSG_IFprotIntegrated_common, "VSG_IFproteinIntegrated_common.csv")
write.csv(VSG_IFprotIntegrated_opposite, "VSG_IFproteinIntegrated_opposite.csv")
write.csv(VSG_IFprotIntegrated, "VSG_IFproteinIntegrated.csv")


###############################################################################
#
#       Integrate IF, TRF and VGS
#
###############################################################################

library(plyr)
TRF_IF_VSGrna <- join_all(list(rnaTRF, rnaIF,rnaVSG),
                          by='GeneName', 
                          type='left'
)
TRF_IF_VSGprot <- join_all(list(protTRF, protIF, protVSG),
                           by='GeneName',
                           type='left'
)
TRF_IF_VSGrna = TRF_IF_VSGrna[complete.cases(TRF_IF_VSGrna), ]
TRF_IF_VSGprot = TRF_IF_VSGprot[complete.cases(TRF_IF_VSGprot), ]

TRF_IF_VSGrna = TRF_IF_VSGrna[!duplicated(TRF_IF_VSGrna[,1]),]
TRF_IF_VSGprot = TRF_IF_VSGprot[!duplicated(TRF_IF_VSGprot[,1]),]
colnames(TRF_IF_VSGrna) = colnames(TRF_IF_VSGprot) = c("GeneName","TRF_FC","TRF_pvalue","IF_FC","IF_pvalue","VSG_FC","VSG_pvalue")

TRF_IF_VSGrna["Log10PvalTRF_signFC"] = -log10(TRF_IF_VSGrna$TRF_pvalue)*sign(TRF_IF_VSGrna$TRF_FC)
TRF_IF_VSGrna["Log10PvalIF_signFC"] = -log10(TRF_IF_VSGrna$IF_pvalue)*sign(TRF_IF_VSGrna$IF_FC)
TRF_IF_VSGrna["Log10PvalVSG_signFC"] = -log10(TRF_IF_VSGrna$VSG_pvalue)*sign(TRF_IF_VSGrna$VSG_FC)

TRF_IF_VSGprot["Log10PvalTRF_signFC"] = -log10(TRF_IF_VSGprot$TRF_pvalue)*sign(TRF_IF_VSGprot$TRF_FC)
TRF_IF_VSGprot["Log10PvalIF_signFC"] = -log10(TRF_IF_VSGprot$IF_pvalue)*sign(TRF_IF_VSGprot$IF_FC)
TRF_IF_VSGprot["Log10PvalVSG_signFC"] = -log10(TRF_IF_VSGprot$VSG_pvalue)*sign(TRF_IF_VSGprot$VSG_FC)

cor(TRF_IF_VSGrna$Log10PvalVSG_signFC,
    TRF_IF_VSGrna$Log10PvalTRF_signFC, 
    method="p"
)
cor(TRF_IF_VSGrna$Log10PvalVSG_signFC,
    TRF_IF_VSGrna$Log10PvalIF_signFC,
    method="p"
)

cor(TRF_IF_VSGprot$Log10PvalVSG_signFC,
    TRF_IF_VSGprot$Log10PvalTRF_signFC, 
    method="p"
)
cor(TRF_IF_VSGprot$Log10PvalVSG_signFC,
    TRF_IF_VSGprot$Log10PvalIF_signFC,
    method="p"
)

#pdf("Katarina.pdf") 
par(mfrow=c(1,2))
plot(TRF_IF_VSGrna$Log10PvalVSG_signFC,
     TRF_IF_VSGrna$Log10PvalTRF_signFC,
     pch=20,
     col="blue",
     main="Transcriptomics",
     ylab="IF and TRF",
     xlab="VSG"
)
points(TRF_IF_VSGrna$Log10PvalVSG_signFC,
       TRF_IF_VSGrna$Log10PvalIF_signFC,
       pch=18,
       col="red"
)
abline(v=c(-log10(cutoff),log10(cutoff)),
       h=c(-log10(cutoff),log10(cutoff)),
       lty=2
)
text(-5.5, 18,
     expression(rho == 0.17),
     cex=1, 
     col="blue"
)
text(-5.5, 19,
     expression(rho == -0.06), 
     cex=1, 
     col="red"
)
legend("topright", 
       legend=c("IF", "TRF"),
       col=c("red", "blue"),
       pch=c(18,20),
       cex=0.8
)

plot(TRF_IF_VSGprot$Log10PvalVSG_signFC,
     TRF_IF_VSGprot$Log10PvalTRF_signFC,
     pch=20,col="blue",
     main="Proteomics",
     ylab="IF and TRF",
     xlab="VSG"
)
points(TRF_IF_VSGprot$Log10PvalVSG_signFC,
       TRF_IF_VSGprot$Log10PvalIF_signFC,
       pch=18,
       col="red"
)
abline(v=c(-log10(cutoff),log10(cutoff)),
       h=c(-log10(cutoff),log10(cutoff)),
       lty=2
)
text(-4.5,
     4,
     expression(rho == 0.31),
     cex=1, 
     col="blue"
)
text(-4.5,
     4.3,
     expression(rho == -0.12),
     cex=1, 
     col="red"
)
legend("topright", 
       legend=c("IF", "TRF"),
       col=c("red", "blue"),
       pch=c(18,20),
       cex=0.8
)
#dev.off() 
###############################################################################
#
#    VSG study - integrated data
#
###############################################################################

####################
#   integrated plot
####################

#VSGfin1=VSGfin1[order(VSGfin1$pvalue.x,VSGfin1$pvalue.y),]
#plot(VSGfin$log2FoldChange.y,
#     VSGfin$log2FoldChange.x, 
#     xlab="Fold change - transciptome",
#     ylab="Fold change - proteome", 
#     main="VSG",
#     pch=19, 
#     xlim=c(-3.5,7.5),
#     col="grey",
#     cex=1.2
#)
#points(VSGfin[VSGfin$pvalue.x < 0.05,]$log2FoldChange.y,
#       VSGfin[VSGfin$pvalue.x < 0.05,]$log2FoldChange.x, 
#       col="darkslategray3",
#       pch=19, 
#       cex=1.2
#)
#points(VSGfin[VSGfin$pvalue.y < 0.05,]$log2FoldChange.y, 
#       VSGfin[VSGfin$pvalue.y < 0.05,]$log2FoldChange.x, 
#       col="plum3", 
#       pch=19,
#       cex=1.2
#)
#points(VSGfin1$log2FoldChange.y,
#       VSGfin1$log2FoldChange.x,
#       col="slateblue3", 
#       pch=19, 
#       cex=1.2
#)
#text(VSGfin1$log2FoldChange.y[1:20],
#     VSGfin1$log2FoldChange.x[1:20], 
#     labels=VSGfin1$GeneName[1:20],
#     cex= 0.9, 
#     pos=3, 
#     col="red",
#     font=2
#)
#text(-3, 
#      1.5,
#      expression(rho == 0.44), 
#      cex=1.5,
#      col="blue"
#)
#abline(v=0,
#       h=0,
#       lty=2
#)

#plot(VSGfin$log10Pval_Trans,
#     VSGfin$log10Pval_Prote, 
#     xlab="-Log10Pval(transcriptome)*sign(FC)",
#     ylab="-Log10Pval(proteome)*sign(FC)", 
#     main="VSG", 
#     pch=19,
#     xlim=c(-6,7.5),
#     col="grey"
#)
#points(VSGfin1$log10Pval_Trans, 
#       VSGfin1$log10Pval_Prote,
#       col="forestgreen", 
#       pch=19
#)
#text(-4.5, 
#     6, 
#     expression(correlaton:rho == 0.44), 
#     cex=0.8,
#     col="blue"
#)
#abline(v=c(-log10(cutoff),log10(cutoff)),
#       h=c(-log10(cutoff),log10(cutoff)),
#       lty=2
#)


VSGfin$significance=ifelse((VSGfin$pvalue.x < 0.05) & (VSGfin$pvalue.y < 0.05), "p-value<0.05", "p-value>0.05")

#setwd("/Users/viktorian.miok/Documents/consultation/Katarina/Integration_Data/Integration/results/VSG/plots/")
#pdf("VSG_integrated.pdf")  
b <- ggplot(VSGfin, aes(x=log10Pval_Trans, y=log10Pval_Prote)) + 
           geom_point(aes(color=significance), show.legend=FALSE) +
           scale_color_manual(values=c("forestgreen", "grey")) +
           geom_hline(yintercept=c(-log10(cutoff),log10(cutoff)),linetype="dashed") +
           geom_vline( xintercept=c(-log10(cutoff),log10(cutoff)),linetype="dashed") +
           xlab("-Log10Pval(transcriptome)*sign(FC)") + 
           ylab("-Log10Pval(proteome)*sign(FC)")+
        annotate(geom="text",
                 x=-5.7,
                 y=6.3, 
                 label=expression(correlaton:rho == 0.44), 
                 color="blue") + 
        ggtitle("VSG study")
#dev.off()
#############################################
#   tables of significant genes and proteins
#############################################
colnames(VSGfin) = colnames(VSGfin1) = c("GeneName","log2FoldChange.protein","pvalue.protein","log2FoldChange.rna",
                                         "pvalue.rna","-Log10Pval(transcriptome)*sign(FC)","-Log10Pval(proteome)*sign(FC)")
VSGintegrated_rna = VSGfin[VSGfin$pvalue.rna < 0.05,]
VSGintegrated_protein = VSGfin[VSGfin$pvalue.protein < 0.05,]
VSGintegrated = VSGfin1
#setwd("/Users/viktorian.miok/Documents/consultation/Katarina/Integration_Data/Integration/results/VSG")
#write.csv(VSGintegrated_rna, "VSG_integrated_RNA.csv")
#write.csv(VSGintegrated_protein, "VSG_integrated_protein.csv")
#write.csv(VSGintegrated, "VSG_integrated_RNA_protein.csv")

########################################
#   enrichment of significant pathways
########################################

sig.gene <- bitr(VSGfin1$GeneName,
                 fromType="SYMBOL",
                 toType="ENTREZID",
                 OrgDb=org.Mm.eg.db
)
kk1 <- enrichKEGG(gene=sig.gene[,2], 
                  organism='mmu'
)
ego1 <- enrichGO(gene=sig.gene[,2], 
                 OrgDb=org.Mm.eg.db, 
                 ont="BP", 
                 pAdjustMethod="BH",
                 readable=TRUE
)
ego2 <- enrichGO(gene=sig.gene[,2], 
                 OrgDb=org.Mm.eg.db, 
                 ont="CC",
                 pAdjustMethod="BH",
                 readable=TRUE
)
ego3 <- enrichGO(gene=sig.gene[,2], 
                 OrgDb=org.Mm.eg.db, 
                 ont="MF",
                 pAdjustMethod="BH",
                 readable=TRUE
)

bp1 = ego1@result
setwd('/Users/viktorian.miok/Documents/consultation/Katarina/Integration_Data/Integration/results/VSG')
###write.csv(bp1, "VSG_GO_BP_integrated.csv")

bp1 = bp1[bp1$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp1[,4]))
bp1 = bp1[(lk > 15) & (lk < 500),]
bp1 = bp1[bp1$Count > 4,]
bp1 = bp1[!bp1$ID %in% noliv[,1],]
bp1 = bp1[!bp1$ID %in% obsolete[,1],]
bp1 = bp1[order(bp1$p.adjust,decreasing=TRUE),]
###write.csv(bp1, "VSG_GO_BP_integrated_a.csv")


rk = character()
for(i in 1:(nrow(bp1)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp1)){
    pk = c(pk,sum(unlist(strsplit(bp1[i,8], "/")) %in% unlist(strsplit(bp1[j,8], "/")))/length(unlist(strsplit(bp1[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp1 = bp1[rk1 == "ostavi",]
##write.csv(bp1, "VSG_GO_BP_integrated_b.csv")

bp2 = ego2@result
#write.csv(bp2, "VSG_GO_CC_integrated.csv")

bp2 = bp2[bp2$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp2[,4]))
bp2 = bp2[(lk > 15) & (lk < 500),]
bp2 = bp2[bp2$Count > 4,]
bp2 = bp2[!bp2$ID %in% noliv[,1],]
bp2 = bp2[!bp2$ID %in% obsolete[,1],]
bp2 = bp2[order(bp2$p.adjust,decreasing=TRUE),]
#write.csv(bp2, "VSG_GO_CC_integrated_a.csv")

rk = character()
for(i in 1:(nrow(bp2)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp2)){
    pk = c(pk,sum(unlist(strsplit(bp2[i,8], "/")) %in% unlist(strsplit(bp2[j,8], "/")))/length(unlist(strsplit(bp2[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp2 = bp2[rk1 == "ostavi",]
#write.csv(bp2, "VSG_GO_CC_integrated_b.csv")


bp3 = ego3@result
#write.csv(bp3, "VSG_GO_MF_integrated.csv")

bp3 = bp3[bp3$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp3[,4]))
bp3 = bp3[(lk > 15) & (lk < 500),]
bp3 = bp3[bp3$Count > 4,]
bp3 = bp3[!bp3$ID %in% noliv[,1],]
bp3 = bp3[!bp3$ID %in% obsolete[,1],]
bp3 = bp3[order(bp3$p.adjust,decreasing=TRUE),]
#write.csv(bp3, "VSG_GO_MF_integrated_a.csv")

rk = character()
for(i in 1:(nrow(bp3)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp3)){
    pk = c(pk,sum(unlist(strsplit(bp3[i,8], "/")) %in% unlist(strsplit(bp3[j,8], "/")))/length(unlist(strsplit(bp3[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp3 = bp3[rk1 == "ostavi",]
#write.csv(bp3, "VSG_GO_MF_integrateed_b.csv")

bp4 = kk1@result
for(i in 1:nrow(bp4)){
  bp4[i,8] <- paste(suppressMessages(bitr(unlist(strsplit(bp4[i,8], "/")),
                                          fromType="ENTREZID",
                                          toType="SYMBOL",
                                          OrgDb=org.Mm.eg.db)[,2]),
                    collapse="/"
)
}

#write.csv(bp4, "VSG_KEGG_intgerated.csv")

bp4 = bp4[bp4$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp4[,4]))
bp4 = bp4[(lk > 15) & (lk < 500),]
bp4 = bp4[bp4$Count > 4,]
bp4 = bp4[!bp4$ID %in% noliv[,1],]
bp4 = bp4[!bp4$ID %in% obsolete[,1],]
bp4 = bp4[order(bp4$p.adjust,decreasing=TRUE),]
#write.csv(bp4, "VSG_KEGG_trimed_integrated_a.csv")

rk = character()
for(i in 1:(nrow(bp4)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp4)){
    pk = c(pk,sum(unlist(strsplit(bp4[i,8], "/")) %in% unlist(strsplit(bp4[j,8], "/")))/length(unlist(strsplit(bp4[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp4 = bp4[rk1 == "ostavi",]
#write.csv(bp4, "VSG_KEGG_integrated_b.csv")

###############################################################################
#
#      TRF study - integrated data
#
###############################################################################

####################
#   integrated plot
####################

#plot(TRFfin$log10Pval_Trans, 
#     TRFfin$log10Pval_Prote,
#     xlab="-Log10Pval(transcriptome)*sign(FC)",
#     ylab="-Log10Pval(proteome)*sign(FC)",
#     main="TRF", 
#     pch=19, 
#     xlim=c(-6,7.5), 
#    col="grey"
#)
#points(TRFfin1$log10Pval_Trans, 
#      TRFfin1$log10Pval_Prote, 
#      col="dodgerblue3", 
#      pch=19
#)
#text(-4, 
#      4,
#      expression(correlaton:rho == 0.24), 
#      cex=0.8,
#      col="blue"
#)
#text(TRFfin1$log10Pval_Trans[TRFfin1$GeneName %in% TRFgen], 
#     TRFfin1$log10Pval_Prote[TRFfin1$GeneName %in% TRFgen],
#     labels=TRFfin1$GeneName[TRFfin1$GeneName %in% TRFgen],
#     cex= 1.2, 
#     pos=3,
#     font=2
#)
#abline(v=c(-log10(cutoff),log10(cutoff)),
#       h=c(-log10(cutoff),log10(cutoff)),
#       lty=2
#)

TRFgen = c("Aadat","Aldh3a2","Aox1","Cat","Ehhadh","Haao","Inmt","Tdo2","Rpl28","Rpl37a","Rpl7a","Rplp0","Rps24","Adh1","Comt","Fah","Hgd")
TRFfin$significance <- ifelse((TRFfin$pvalue.x < 0.05) & (TRFfin$pvalue.y < 0.05),
                              "p-value<0.05",
                              "p-value>0.05"
)
TRFfin$label <- ifelse(TRFfin$GeneName %in% TRFgen,
                       "name",
                       ""
)
for(i in 1:nrow(TRFfin)){
  if(TRFfin$label[i] == "name") TRFfin$label[i] = TRFfin$GeneName[i]
}

#setwd("/Users/viktorian.miok/Documents/consultation/Katarina/Integration_Data/Integration/results/TRF/plots/")
#pdf("TRF_integrated.pdf") 
b <- ggplot(TRFfin, aes(x=log10Pval_Trans, y=log10Pval_Prote)) + 
           geom_point(aes(color=significance), show.legend=FALSE) +
           scale_color_manual(values=c("dodgerblue3", "grey")) +
           geom_hline(yintercept=c(-log10(cutoff),log10(cutoff)),linetype="dashed") +
           geom_vline( xintercept=c(-log10(cutoff),log10(cutoff)),linetype="dashed") +
           xlab("-Log10Pval(transcriptome)*sign(FC)") +
           ylab("-Log10Pval(proteome)*sign(FC)") +
           annotate(geom="text", x=-5, y=5.1, label=expression(correlaton:rho == 0.24), color="blue") +
           ggtitle("TRF study") +
           geom_text_repel(aes(label=label), 
                           size=5,
                           box.padding=unit(2.5, "lines"), 
                           point.padding=unit(0.1, "lines"),
                           segment.size=0.15
           ) 
#dev.off()

#############################################
#   tables of significant genes and proteins
#############################################
colnames(TRFfin)=colnames(TRFfin1) = c("GeneName","log2FoldChange.protein","pvalue.protein","log2FoldChange.rna",
                                       "pvalue.rna","-Log10Pval(transcriptome)*sign(FC)","-Log10Pval(proteome)*sign(FC)")
TRFintegrated_rna = TRFfin[TRFfin$pvalue.rna < 0.05,]
TRFintegrated_protein = TRFfin[TRFfin$pvalue.protein < 0.05,]
TRFintegrated = TRFfin1
#setwd("/Users/viktorian.miok/Documents/consultation/Katarina/Integration_Data/Integration/results/TRF")
#write.csv(TRFintegrated_rna, "TRF_integrated_RNA.csv")
#write.csv(TRFintegrated_protein, "TRF_integrated_protein.csv")
#write.csv(TRFintegrated, "TRF_integrated_RNA_protein.csv")

########################################
#   enrichment of significant pathways
########################################

sig.gene <- bitr(TRFfin1$GeneName,
                 fromType="SYMBOL",
                 toType="ENTREZID",
                 OrgDb=org.Mm.eg.db
)
kk1 <- enrichKEGG(gene=sig.gene[,2], 
                  organism='mmu'
)
ego1 <- enrichGO(gene=sig.gene[,2], 
                 OrgDb=org.Mm.eg.db,
                 ont="BP", 
                 pAdjustMethod="BH",
                 readable=TRUE
)
ego2 <- enrichGO(gene=sig.gene[,2],
                 OrgDb=org.Mm.eg.db, 
                 ont="CC",
                 pAdjustMethod="BH",
                 readable=TRUE
)
ego3 <- enrichGO(gene=sig.gene[,2], 
                 OrgDb=org.Mm.eg.db, 
                 ont="MF", 
                 pAdjustMethod="BH",
                 readable=TRUE
)

bp1 = ego1@result
setwd('/Users/viktorian.miok/Documents/consultation/Katarina/Integration_Data/Integration/results/TRF')
#write.csv(bp1, "TRF_GO_BP_integrated.csv")

bp1 = bp1[bp1$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp1[,4]))
bp1 = bp1[(lk > 15) & (lk < 500),]
bp1 = bp1[bp1$Count > 4,]
bp1 = bp1[!bp1$ID %in% noliv[,1],]
bp1 = bp1[!bp1$ID %in% obsolete[,1],]
bp1 = bp1[order(bp1$p.adjust,decreasing=TRUE),]
#write.csv(bp1, "TRF_GO_BP_integrated_a.csv")

rk = character()
for(i in 1:(nrow(bp1)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp1)){
    pk = c(pk,sum(unlist(strsplit(bp1[i,8], "/")) %in% unlist(strsplit(bp1[j,8], "/")))/length(unlist(strsplit(bp1[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp1 = bp1[rk1 == "ostavi",]
#write.csv(bp1, "TRF_GO_BP_integrated_b.csv")

bp2 = ego2@result
#write.csv(bp2, "TRF_GO_CC_integrated.csv")

bp2 = bp2[bp2$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp2[,4]))
bp2 = bp2[(lk > 15) & (lk < 500),]
bp2 = bp2[bp2$Count > 4,]
bp2 = bp2[!bp2$ID %in% noliv[,1],]
bp2 = bp2[!bp2$ID %in% obsolete[,1],]
bp2 = bp2[order(bp2$p.adjust,decreasing=TRUE),]
#write.csv(bp2, "TRF_GO_CC_integrated_a.csv")

rk = character()
for(i in 1:(nrow(bp2)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp2)){
    pk = c(pk,sum(unlist(strsplit(bp2[i,8], "/")) %in% unlist(strsplit(bp2[j,8], "/")))/length(unlist(strsplit(bp2[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp2 = bp2[rk1 == "ostavi",]
#write.csv(bp2, "TRF_GO_CC_integrated_b.csv")

bp3 = ego3@result
#write.csv(bp3, "TRF_GO_MF_integrated.csv")

bp3 = bp3[bp3$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp3[,4]))
bp3 = bp3[(lk > 15) & (lk < 500),]
bp3 = bp3[bp3$Count > 4,]
bp3 = bp3[!bp3$ID %in% noliv[,1],]
bp3 = bp3[!bp3$ID %in% obsolete[,1],]
bp3 = bp3[order(bp3$p.adjust,decreasing=TRUE),]
#write.csv(bp3, "TRF_GO_MF_integrated_a.csv")

rk = character()
for(i in 1:(nrow(bp3)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp3)){
    pk = c(pk,sum(unlist(strsplit(bp3[i,8], "/")) %in% unlist(strsplit(bp3[j,8], "/")))/length(unlist(strsplit(bp3[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp3 = bp3[rk1 == "ostavi",]
#write.csv(bp3, "TRF_GO_MF_integrated_b.csv")

bp4 = kk1@result
for(i in 1:nrow(bp4)){
  bp4[i,8] <- paste(suppressMessages(bitr(unlist(strsplit(bp4[i,8], "/")),
                                          fromType="ENTREZID",
                                          toType="SYMBOL",
                                          OrgDb=org.Mm.eg.db)[,2]),
                    collapse="/"
)
}

#write.csv(bp4, "TRF_KEGG_integrated.csv")

bp4 = bp4[bp4$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp4[,4]))
bp4 = bp4[(lk > 15) & (lk < 500),]
bp4 = bp4[bp4$Count > 4,]
bp4 = bp4[!bp4$ID %in% noliv[,1],]
bp4 = bp4[!bp4$ID %in% obsolete[,1],]
bp4 = bp4[order(bp4$p.adjust,decreasing=TRUE),]
#write.csv(bp4, "TRF_KEGG_integrated_a.csv")

rk = character()
for(i in 1:(nrow(bp4)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp4)){
    pk = c(pk,sum(unlist(strsplit(bp4[i,8], "/")) %in% unlist(strsplit(bp4[j,8], "/")))/length(unlist(strsplit(bp4[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp4 = bp4[rk1 == "ostavi",]
#write.csv(bp4, "TRF_KEGG_integrated_b.csv")

###############################################################################
#
#     IF study - integrated data
#
###############################################################################

#plot(IFfin$log10Pval_Trans, 
#     IFfin$log10Pval_Prote, 
#     xlab="-Log10Pval(transcriptome)*sign(FC)",
#     ylab="-Log10Pval(proteome)*sign(FC)",
#     main="IF",
#     pch=19,
#     xlim=c(-6,7.5),
#     col="grey"
#)
#points(IFfin1$log10Pval_Trans,
#       IFfin1$log10Pval_Prote,
#       col="firebrick", 
#       pch=19)
#text(-4,
#     4.5, 
#     expression(correlaton:rho == 0.19), 
#     cex=0.8, 
#     col="blue"
#)
#text(IFfin1$log10Pval_Trans[IFfin1$GeneName %in% IFgen],
#     IFfin1$log10Pval_Prote[IFfin1$GeneName %in% IFgen], 
#     labels=IFfin1$GeneName[IFfin1$GeneName %in% IFgen],
#     cex= 0.9, 
#     pos=3, 
#     col="red", 
#     font=2
#)
#abline(v=c(-log10(cutoff),log10(cutoff)),h=c(-log10(cutoff),log10(cutoff)),lty=2)

IFgen = c("Baat","Fads1","Fads2","Scd1","Crot","Hacl1","Mlycd","Pmvk","Cyc1","Ndufa3","Ndufa9")
IFfin$significance <- ifelse((IFfin$pvalue.x < 0.05) & (IFfin$pvalue.y < 0.05),
                             "p-value<0.05",
                             "p-value>0.05"
)
IFfin$label <- ifelse(IFfin$GeneName %in% IFgen,
                      "name",
                      ""
)
for(i in 1:nrow(IFfin)){
  if(IFfin$label[i] == "name") IFfin$label[i] = IFfin$GeneName[i]
}

setwd("/Users/viktorian.miok/Documents/consultation/Katarina/Integration_Data/Integration/results/IF/plots/")
pdf("IF_integrated.pdf") 
b <- ggplot(IFfin, aes(x=log10Pval_Trans, y=log10Pval_Prote)) +
           geom_point(aes(color=significance), show.legend=FALSE) +
           scale_color_manual(values=c("firebrick", "grey")) +
           geom_hline(yintercept=c(-log10(cutoff), log10(cutoff)), linetype="dashed") +
           geom_vline( xintercept=c(-log10(cutoff), log10(cutoff)), linetype="dashed") +
           xlab("-Log10Pval(transcriptome)*sign(FC)") + 
           ylab("-Log10Pval(proteome)*sign(FC)")+
           annotate(geom="text",
                    x=-7.5,
                    y=6,
                    label=expression(correlaton:rho == 0.19),
                    color="blue") + 
           ggtitle("TRF study")+
           geom_text_repel(aes(label=label), 
                           size=5,
                           box.padding=unit(2, "lines"), 
                           point.padding=unit(0.1, "lines"),
                           segment.size=0.15
           ) 
dev.off()

#############################################
#   tables of significant genes and proteins
#############################################
colnames(IFfin) = colnames(IFfin1) = c("GeneName","log2FoldChange.protein","pvalue.protein","log2FoldChange.rna",
                                       "pvalue.rna","-Log10Pval(transcriptome)*sign(FC)","-Log10Pval(proteome)*sign(FC)")
IFintegrated_rna = IFfin[IFfin$pvalue.rna < 0.05,]
IFintegrated_protein = IFfin[IFfin$pvalue.protein < 0.05,]
IFintegrated = IFfin1
#setwd("/Users/viktorian.miok/Documents/consultation/Katarina/Integration_Data/Integration/results/IF")
#write.csv(IFintegrated_rna, "IF_integrated_RNA.csv")
#write.csv(IFintegrated_protein, "IF_integrated_protein.csv")
#write.csv(IFintegrated, "IF_integrated_RNA_protein.csv")

########################################
#   enrichment of significant pathways
########################################

sig.gene <- bitr(IFfin1$GeneName,
                 fromType="SYMBOL",
                 toType="ENTREZID",
                 OrgDb=org.Mm.eg.db
)
kk1 <- enrichKEGG(gene=sig.gene[,2], 
                  organism='mmu'
)
ego1 <- enrichGO(gene=sig.gene[,2], 
                 OrgDb=org.Mm.eg.db,
                 ont="BP",
                 pAdjustMethod="BH",
                 readable=TRUE
)
ego2 <- enrichGO(gene=sig.gene[,2], 
                 OrgDb=org.Mm.eg.db, 
                 ont="CC",
                 pAdjustMethod="BH",
                 readable=TRUE
)
ego3 <- enrichGO(gene=sig.gene[,2],
                 OrgDb=org.Mm.eg.db,
                 ont="MF",
                 pAdjustMethod="BH", 
                 readable=TRUE
)

bp1 = ego1@result
setwd('/Users/viktorian.miok/Documents/consultation/Katarina/Integration_Data/Integration/results/IF')
#write.csv(bp1, "IF_GO_BP_integrated.csv")

bp1 = bp1[bp1$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp1[,4]))
bp1 = bp1[(lk > 15) & (lk < 500),]
bp1 = bp1[bp1$Count > 4,]
bp1 = bp1[!bp1$ID %in% noliv[,1],]
bp1 = bp1[!bp1$ID %in% obsolete[,1],]
bp1 = bp1[order(bp1$p.adjust,decreasing=TRUE),]
#write.csv(bp1, "IF_GO_BP_integrated_a.csv")

rk = character()
for(i in 1:(nrow(bp1)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp1)){
    pk = c(pk,sum(unlist(strsplit(bp1[i,8], "/")) %in% unlist(strsplit(bp1[j,8], "/")))/length(unlist(strsplit(bp1[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp1 = bp1[rk1 == "ostavi",]
#write.csv(bp1, "IF_GO_BP_integrated_b.csv")

bp2 = ego2@result
#write.csv(bp2, "IF_GO_CC_integrated.csv")

bp2 = bp2[bp2$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp2[,4]))
bp2 = bp2[(lk > 15) & (lk < 500),]
bp2 = bp2[bp2$Count > 4,]
bp2 = bp2[!bp2$ID %in% noliv[,1],]
bp2 = bp2[!bp2$ID %in% obsolete[,1],]
bp2 = bp2[order(bp2$p.adjust,decreasing=TRUE),]
#write.csv(bp2, "IF_GO_CC_integrated_a.csv")

rk = character()
for(i in 1:(nrow(bp2)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp2)){
    pk = c(pk,sum(unlist(strsplit(bp2[i,8], "/")) %in% unlist(strsplit(bp2[j,8], "/")))/length(unlist(strsplit(bp2[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp2 = bp2[rk1 == "ostavi",]
#write.csv(bp2, "IF_GO_CC_integrated_b.csv")

bp3 = ego3@result
#write.csv(bp3, "IF_GO_MF_integrated.csv")

bp3 = bp3[bp3$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp3[,4]))
bp3 = bp3[(lk > 15) & (lk < 500),]
bp3 = bp3[bp3$Count > 4,]
bp3 = bp3[!bp3$ID %in% noliv[,1],]
bp3 = bp3[!bp3$ID %in% obsolete[,1],]
bp3 = bp3[order(bp3$p.adjust,decreasing=TRUE),]
#write.csv(bp3, "IF_GO_MF_integrated_a.csv")

rk = character()
for(i in 1:(nrow(bp3)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp3)){
    pk = c(pk,sum(unlist(strsplit(bp3[i,8], "/")) %in% unlist(strsplit(bp3[j,8], "/")))/length(unlist(strsplit(bp3[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp3 = bp3[rk1 == "ostavi",]
#write.csv(bp3, "IF_GO_MF_integrated_b.csv")

bp4 = kk1@result
for(i in 1:nrow(bp4)){
  bp4[i,8] <- paste(suppressMessages(bitr(unlist(strsplit(bp4[i,8], "/")),
                                          fromType="ENTREZID",
                                          toType="SYMBOL",
                                          OrgDb=org.Mm.eg.db)[,2]),
                    collapse="/"
)
}

#write.csv(bp4, "IF_KEGG_integrated.csv")

bp4 = bp4[bp4$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp4[,4]))
bp4 = bp4[(lk > 15) & (lk < 500),]
bp4 = bp4[bp4$Count > 4,]
bp4 = bp4[!bp4$ID %in% noliv[,1],]
bp4 = bp4[!bp4$ID %in% obsolete[,1],]
bp4 = bp4[order(bp4$p.adjust,decreasing=TRUE),]
#write.csv(bp4, "IF_KEGG_integrated_a.csv")

rk = character()
for(i in 1:(nrow(bp4)- 1)){
  pk = numeric()
  for(j in (i+1):nrow(bp4)){
    pk = c(pk,sum(unlist(strsplit(bp4[i,8], "/")) %in% unlist(strsplit(bp4[j,8], "/")))/length(unlist(strsplit(bp4[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp4 = bp4[rk1 == "ostavi",]
#write.csv(bp4, "IF_KEGG_integrated_b.csv")

###############################################################################
#
#    Poroteomics data enrichment
#
###############################################################################

setwd('/Users/viktorian.miok/Documents/consultation/Katarina/Integration_Data/RNA_seq_data')
load("ProteomOutput.RData")
###############################################################################
# IF - up
###############################################################################
setwd('/Users/viktorian.miok/Documents/consultation/Katarina/Integration_Data/Proteomics_data/results/IF')

sig.gene <- bitr(ProteomOutput$diffIF[(ProteomOutput$diffIF$log2FoldChange > 0) & (ProteomOutput$diffIF$pvalue < 0.05),1],
                 fromType="SYMBOL",
                 toType="ENTREZID",
                 OrgDb=org.Mm.eg.db
)
kk1 <- enrichKEGG(gene=sig.gene[,2],
                  organism='mmu'
)
ego1 <- enrichGO(gene=sig.gene[,2],
                 OrgDb=org.Mm.eg.db,
                 ont="BP", 
                 pAdjustMethod="BH", 
                 readable=TRUE
)
ego2 <- enrichGO(gene=sig.gene[,2],
                 OrgDb=org.Mm.eg.db,
                 ont="CC",
                 pAdjustMethod="BH", 
                 readable=TRUE
)
ego3 <- enrichGO(gene=sig.gene[,2], 
                 OrgDb=org.Mm.eg.db, 
                 ont="MF", 
                 pAdjustMethod="BH",
                 readable=TRUE
)

bp1 = ego1@result
#write.csv(bp1, "IF_GO_BP_up_prot.csv")

bp1 = bp1[bp1$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp1[,4]))
bp1 = bp1[(lk > 15) & (lk < 500),]
bp1 = bp1[bp1$Count > 4,]
bp1 = bp1[!bp1$ID %in% noliv[,1],]
bp1 = bp1[!bp1$ID %in% obsolete[,1],]
bp1 = bp1[order(bp1$p.adjust,decreasing=TRUE),]
#write.csv(bp1, "IF_GO_BP_trimed_up_prot_a.csv")

rk = character()
for(i in 1:(nrow(bp1)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp1)){
    pk = c(pk,sum(unlist(strsplit(bp1[i,8], "/")) %in% unlist(strsplit(bp1[j,8], "/")))/length(unlist(strsplit(bp1[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp1 = bp1[rk1 == "ostavi",]
#write.csv(bp1, "IF_GO_BP_trimed_up_prot_b.csv")

bp2 = ego2@result
#write.csv(bp2, "IF_GO_CC_up_prot.csv")

bp2 = bp2[bp2$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp2[,4]))
bp2 = bp2[(lk > 15) & (lk < 500),]
bp2 = bp2[bp2$Count > 4,]
bp2 = bp2[!bp2$ID %in% noliv[,1],]
bp2 = bp2[!bp2$ID %in% obsolete[,1],]
bp2 = bp2[order(bp2$p.adjust,decreasing=TRUE),]
#write.csv(bp2, "IF_GO_CC_trimed_up_prot_a.csv")

rk = character()
for(i in 1:(nrow(bp2)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp2)){
    pk = c(pk,sum(unlist(strsplit(bp2[i,8], "/")) %in% unlist(strsplit(bp2[j,8], "/")))/length(unlist(strsplit(bp2[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp2 = bp2[rk1 == "ostavi",]
#write.csv(bp2, "IF_GO_CC_trimed_up_prot_b.csv")

bp3 = ego3@result
#write.csv(bp3, "IF_GO_MF_up_prot.csv")

bp3 = bp3[bp3$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp3[,4]))
bp3 = bp3[(lk > 15) & (lk < 500),]
bp3 = bp3[bp3$Count > 4,]
bp3 = bp3[!bp3$ID %in% noliv[,1],]
bp3 = bp3[!bp3$ID %in% obsolete[,1],]
bp3 = bp3[order(bp3$p.adjust,decreasing=TRUE),]
#write.csv(bp3, "IF_GO_MF_trimed_up_prot_a.csv")

rk = character()
for(i in 1:(nrow(bp3)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp3)){
    pk = c(pk,sum(unlist(strsplit(bp3[i,8], "/")) %in% unlist(strsplit(bp3[j,8], "/")))/length(unlist(strsplit(bp3[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp3 = bp3[rk1 == "ostavi",]
#write.csv(bp3, "IF_GO_MF_trimed_up_prot_b.csv")

bp4 = kk1@result
for(i in 1:nrow(bp4)){
  bp4[i,8] <- paste(suppressMessages(bitr(unlist(strsplit(bp4[i,8], "/")),
                                          fromType="ENTREZID",
                                          toType="SYMBOL",
                                          OrgDb=org.Mm.eg.db)[,2]),
                    collapse="/"
)
}

#write.csv(bp4, "IF_KEGG_up_prot.csv")

bp4 = bp4[bp4$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp4[,4]))
bp4 = bp4[(lk > 15) & (lk < 500),]
bp4 = bp4[bp4$Count > 4,]
bp4 = bp4[!bp4$ID %in% noliv[,1],]
bp4 = bp4[!bp4$ID %in% obsolete[,1],]
bp4 = bp4[order(bp4$p.adjust,decreasing=TRUE),]

#write.csv(bp4, "IF_KEGG_trimed_up_prot_a.csv")

rk = character()
for(i in 1:(nrow(bp4)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp4)){
    pk = c(pk,sum(unlist(strsplit(bp4[i,8], "/")) %in% unlist(strsplit(bp4[j,8], "/")))/length(unlist(strsplit(bp4[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp4 = bp4[rk1 == "ostavi",]
#write.csv(bp4, "IF_KEGG_trimed_up_prot_b.csv")


###############################################################################
# IF - down
###############################################################################
sig.gene <- bitr(ProteomOutput$diffIF[(ProteomOutput$diffIF$log2FoldChange < 0) & (ProteomOutput$diffIF$pvalue < 0.05),1],
                 fromType="SYMBOL",
                 toType="ENTREZID",
                 OrgDb=org.Mm.eg.db
)
kk1 <- enrichKEGG(gene=sig.gene[,2], 
                  organism='mmu'
)
ego1 <- enrichGO(gene=sig.gene[,2],
                 OrgDb=org.Mm.eg.db,
                 ont="BP",
                 pAdjustMethod="BH",
                 readable=TRUE
)
ego2 <- enrichGO(gene=sig.gene[,2],
                 OrgDb=org.Mm.eg.db,
                 ont="CC", 
                 pAdjustMethod="BH", 
                 readable=TRUE
)
ego3 <- enrichGO(gene=sig.gene[,2], 
                 OrgDb=org.Mm.eg.db,
                 ont="MF", 
                 pAdjustMethod="BH", 
                 readable=TRUE
)
bp1 = ego1@result
#write.csv(bp1, "IF_GO_BP_down_prot.csv")

bp1 = bp1[bp1$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp1[,4]))
bp1 = bp1[(lk > 15) & (lk < 500),]
bp1 = bp1[bp1$Count > 4,]
bp1 = bp1[!bp1$ID %in% noliv[,1],]
bp1 = bp1[!bp1$ID %in% obsolete[,1],]
bp1 = bp1[order(bp1$p.adjust,decreasing=TRUE),]
#write.csv(bp1, "IF_GO_BP_trimed_down_prot_a.csv")

rk = character()
for(i in 1:(nrow(bp1)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp1)){
    pk = c(pk,sum(unlist(strsplit(bp1[i,8], "/")) %in% unlist(strsplit(bp1[j,8], "/")))/length(unlist(strsplit(bp1[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp1 = bp1[rk1 == "ostavi",]
#write.csv(bp1, "IF_GO_BP_trimed_down_prot_b.csv")

bp2 = ego2@result
#write.csv(bp2, "IF_GO_CC_down_prot.csv")

bp2 = bp2[bp2$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp2[,4]))
bp2 = bp2[(lk > 15) & (lk < 500),]
bp2 = bp2[bp2$Count > 4,]
bp2 = bp2[!bp2$ID %in% noliv[,1],]
bp2 = bp2[!bp2$ID %in% obsolete[,1],]
bp2 = bp2[order(bp2$p.adjust,decreasing=TRUE),]
#write.csv(bp2, "IF_GO_CC_trimed_down_prot_a.csv")

rk = character()
for(i in 1:(nrow(bp2)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp2)){
    pk = c(pk,sum(unlist(strsplit(bp2[i,8], "/")) %in% unlist(strsplit(bp2[j,8], "/")))/length(unlist(strsplit(bp2[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp2 = bp2[rk1 == "ostavi",]
#write.csv(bp2, "IF_GO_CC_trimed_down_prot_b.csv")

bp3 = ego3@result
#write.csv(bp3, "IF_GO_MF_down_prot.csv")

bp3 = bp3[bp3$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp3[,4]))
bp3 = bp3[(lk > 15) & (lk < 500),]
bp3 = bp3[bp3$Count > 4,]
bp3 = bp3[!bp3$ID %in% noliv[,1],]
bp3 = bp3[!bp3$ID %in% obsolete[,1],]
bp3 = bp3[order(bp3$p.adjust,decreasing=TRUE),]
#write.csv(bp3, "IF_GO_MF_trimed_down_prot_a.csv")

rk = character()
for(i in 1:(nrow(bp3)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp3)){
    pk = c(pk,sum(unlist(strsplit(bp3[i,8], "/")) %in% unlist(strsplit(bp3[j,8], "/")))/length(unlist(strsplit(bp3[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp3 = bp3[rk1 == "ostavi",]
#write.csv(bp3, "IF_GO_MF_trimed_down_prot_b.csv")

bp4 = kk1@result
for(i in 1:nrow(bp4)){
  bp4[i,8] <- paste(suppressMessages(bitr(unlist(strsplit(bp4[i,8], "/")),
                                          fromType="ENTREZID",
                                          toType="SYMBOL",
                                          OrgDb=org.Mm.eg.db)[,2]),
                    collapse="/"
)
}

#write.csv(bp4, "IF_KEGG_down_prot.csv")

bp4 = bp4[bp4$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp4[,4]))
bp4 = bp4[(lk > 15) & (lk < 500),]
bp4 = bp4[bp4$Count > 4,]
bp4 = bp4[!bp4$ID %in% noliv[,1],]
bp4 = bp4[!bp4$ID %in% obsolete[,1],]
bp4 = bp4[order(bp4$p.adjust,decreasing=TRUE),]

#write.csv(bp4, "IF_KEGG_trimed_down_prot_a.csv")

rk = character()
for(i in 1:(nrow(bp4)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp4)){
    pk = c(pk,sum(unlist(strsplit(bp4[i,8], "/")) %in% unlist(strsplit(bp4[j,8], "/")))/length(unlist(strsplit(bp4[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp4 = bp4[rk1 == "ostavi",]
#write.csv(bp4, "IF_KEGG_trimed_down_prot_b.csv")

###############################################################################
# VSG - up
###############################################################################
setwd('/Users/viktorian.miok/Documents/consultation/Katarina/Integration_Data/Proteomics_data/results/VSG')

sig.gene <- bitr(ProteomOutput$diffVSG[(ProteomOutput$diffVSG$log2FoldChange > 0) & (ProteomOutput$diffVSG$pvalue < 0.05),1],
                 fromType="SYMBOL",
                 toType="ENTREZID",
                 OrgDb=org.Mm.eg.db
)
kk1 <- enrichKEGG(gene=sig.gene[,2], 
                  organism='mmu'
)
ego1 <- enrichGO(gene=sig.gene[,2],
                 OrgDb=org.Mm.eg.db, 
                 ont="BP",
                 pAdjustMethod="BH",
                 readable=TRUE
)
ego2 <- enrichGO(gene=sig.gene[,2],
                 OrgDb=org.Mm.eg.db,
                 ont="CC", 
                 pAdjustMethod="BH", 
                 readable=TRUE
)
ego3 <- enrichGO(gene=sig.gene[,2],
                 OrgDb=org.Mm.eg.db,
                 ont="MF",
                 pAdjustMethod="BH",
                 readable=TRUE)

bp1 = ego1@result
#write.csv(bp1, "VSG_GO_BP_up_prot.csv")

bp1 = bp1[bp1$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp1[,4]))
bp1 = bp1[(lk > 15) & (lk < 500),]
bp1 = bp1[bp1$Count > 4,]
bp1 = bp1[!bp1$ID %in% noliv[,1],]
bp1 = bp1[!bp1$ID %in% obsolete[,1],]
bp1 = bp1[order(bp1$p.adjust,decreasing=TRUE),]
#write.csv(bp1, "VSG_GO_BP_trimed_up_prot_a.csv")

rk = character()
for(i in 1:(nrow(bp1)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp1)){
    pk = c(pk,sum(unlist(strsplit(bp1[i,8], "/")) %in% unlist(strsplit(bp1[j,8], "/")))/length(unlist(strsplit(bp1[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp1 = bp1[rk1 == "ostavi",]
#write.csv(bp1, "VSG_GO_BP_trimed_up_prot_b.csv")

bp2 = ego2@result
#write.csv(bp2, "VSG_GO_CC_up_prot.csv")

bp2 = bp2[bp2$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp2[,4]))
bp2 = bp2[(lk > 15) & (lk < 500),]
bp2 = bp2[bp2$Count > 4,]
bp2 = bp2[!bp2$ID %in% noliv[,1],]
bp2 = bp2[!bp2$ID %in% obsolete[,1],]
bp2 = bp2[order(bp2$p.adjust,decreasing=TRUE),]
#write.csv(bp2, "VSG_GO_CC_trimed_up_prot_a.csv")

rk = character()
for(i in 1:(nrow(bp2)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp2)){
    pk = c(pk,sum(unlist(strsplit(bp2[i,8], "/")) %in% unlist(strsplit(bp2[j,8], "/")))/length(unlist(strsplit(bp2[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp2 = bp2[rk1 == "ostavi",]
#write.csv(bp2, "VSG_GO_CC_trimed_up_prot_b.csv")

bp3 = ego3@result
#write.csv(bp3, "VSG_GO_MF_up_prot.csv")

bp3 = bp3[bp3$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp3[,4]))
bp3 = bp3[(lk > 15) & (lk < 500),]
bp3 = bp3[bp3$Count > 4,]
bp3 = bp3[!bp3$ID %in% noliv[,1],]
bp3 = bp3[!bp3$ID %in% obsolete[,1],]
bp3 = bp3[order(bp3$p.adjust,decreasing=TRUE),]
#write.csv(bp3, "VSG_GO_MF_trimed_up_prot_a.csv")

rk = character()
for(i in 1:(nrow(bp3)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp3)){
    pk = c(pk,sum(unlist(strsplit(bp3[i,8], "/")) %in% unlist(strsplit(bp3[j,8], "/")))/length(unlist(strsplit(bp3[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp3 = bp3[rk1 == "ostavi",]
#write.csv(bp3, "VSG_GO_MF_trimed_up_prot_b.csv")

bp4 = kk1@result
for(i in 1:nrow(bp4)){
  bp4[i,8] <- paste(suppressMessages(bitr(unlist(strsplit(bp4[i,8], "/")),
                                          fromType="ENTREZID",
                                          toType="SYMBOL",
                                          OrgDb=org.Mm.eg.db)[,2]),
                    collapse="/"
  )
}

#write.csv(bp4, "VSG_KEGG_up_prot.csv")

bp4 = bp4[bp4$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp4[,4]))
bp4 = bp4[(lk > 15) & (lk < 500),]
bp4 = bp4[bp4$Count > 4,]
bp4 = bp4[!bp4$ID %in% noliv[,1],]
bp4 = bp4[!bp4$ID %in% obsolete[,1],]
bp4 = bp4[order(bp4$p.adjust,decreasing=TRUE),]

#write.csv(bp4, "VSG_KEGG_trimed_up_prot_a.csv")

rk = character()
for(i in 1:(nrow(bp4)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp4)){
    pk = c(pk,sum(unlist(strsplit(bp4[i,8], "/")) %in% unlist(strsplit(bp4[j,8], "/")))/length(unlist(strsplit(bp4[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp4 = bp4[rk1 == "ostavi",]
#write.csv(bp4, "VSG_KEGG_trimed_up_prot_b.csv")


###############################################################################
# VSG - down
###############################################################################
sig.gene <- bitr(ProteomOutput$diffVSG[(ProteomOutput$diffVSG$log2FoldChange < 0) & (ProteomOutput$diffVSG$pvalue < 0.05), 1],
                 fromType="SYMBOL",
                 toType="ENTREZID",
                 OrgDb=org.Mm.eg.db
)
kk1 <- enrichKEGG(gene=sig.gene[,2], 
                  organism='mmu'
)
ego1 <- enrichGO(gene=sig.gene[,2],
                 OrgDb=org.Mm.eg.db,
                 ont="BP",
                 pAdjustMethod="BH", 
                 readable=TRUE
)
ego2 <- enrichGO(gene=sig.gene[,2],
                 OrgDb=org.Mm.eg.db, 
                 ont="CC",
                 pAdjustMethod="BH",
                 readable=TRUE
)
ego3 <- enrichGO(gene=sig.gene[,2],
                 OrgDb=org.Mm.eg.db,
                 ont="MF",
                 pAdjustMethod="BH",
                 readable=TRUE
)
bp1 = ego1@result
#write.csv(bp1, "VSG_GO_BP_down_prot.csv")

bp1 = bp1[bp1$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp1[,4]))
bp1 = bp1[(lk > 15) & (lk < 500),]
bp1 = bp1[bp1$Count > 4,]
bp1 = bp1[!bp1$ID %in% noliv[,1],]
bp1 = bp1[!bp1$ID %in% obsolete[,1],]
bp1 = bp1[order(bp1$p.adjust,decreasing=TRUE),]
#write.csv(bp1, "VSG_GO_BP_trimed_down_prot_a.csv")

rk = character()
for(i in 1:(nrow(bp1)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp1)){
    pk = c(pk,sum(unlist(strsplit(bp1[i,8], "/")) %in% unlist(strsplit(bp1[j,8], "/")))/length(unlist(strsplit(bp1[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp1 = bp1[rk1 == "ostavi",]
#write.csv(bp1, "VSG_GO_BP_trimed_down_prot_b.csv")

bp2 = ego2@result
#write.csv(bp2, "VSG_GO_CC_down_prot.csv")

bp2 = bp2[bp2$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp2[,4]))
bp2 = bp2[(lk > 15) & (lk < 500),]
bp2 = bp2[bp2$Count > 4,]
bp2 = bp2[!bp2$ID %in% noliv[,1],]
bp2 = bp2[!bp2$ID %in% obsolete[,1],]
bp2 = bp2[order(bp2$p.adjust,decreasing=TRUE),]
#write.csv(bp2, "VSG_GO_CC_trimed_down_prot_a.csv")

rk = character()
for(i in 1:(nrow(bp2)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp2)){
    pk = c(pk,sum(unlist(strsplit(bp2[i,8], "/")) %in% unlist(strsplit(bp2[j,8], "/")))/length(unlist(strsplit(bp2[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp2 = bp2[rk1 == "ostavi",]
#write.csv(bp2, "VSG_GO_CC_trimed_down_prot_b.csv")

bp3 = ego3@result
#write.csv(bp3, "VSG_GO_MF_down_prot.csv")

bp3 = bp3[bp3$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp3[,4]))
bp3 = bp3[(lk > 15) & (lk < 500),]
bp3 = bp3[bp3$Count > 4,]
bp3 = bp3[!bp3$ID %in% noliv[,1],]
bp3 = bp3[!bp3$ID %in% obsolete[,1],]
bp3 = bp3[order(bp3$p.adjust,decreasing=TRUE),]
#write.csv(bp3, "VSG_GO_MF_trimed_down_prot_a.csv")

rk = character()
for(i in 1:(nrow(bp3)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp3)){
    pk = c(pk,sum(unlist(strsplit(bp3[i,8], "/")) %in% unlist(strsplit(bp3[j,8], "/")))/length(unlist(strsplit(bp3[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp3 = bp3[rk1 == "ostavi",]
#write.csv(bp3, "VSG_GO_MF_trimed_down_prot_b.csv")

bp4 = kk1@result
for(i in 1:nrow(bp4)){
  bp4[i,8] <- paste(suppressMessages(bitr(unlist(strsplit(bp4[i,8], "/")),
                                          fromType="ENTREZID",
                                          toType="SYMBOL",
                                          OrgDb=org.Mm.eg.db)[,2]),
                    collapse="/"
)
}

#write.csv(bp4, "VSG_KEGG_down_prot.csv")

bp4 = bp4[bp4$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp4[,4]))
bp4 = bp4[(lk > 15) & (lk < 500),]
bp4 = bp4[bp4$Count > 4,]
bp4 = bp4[!bp4$ID %in% noliv[,1],]
bp4 = bp4[!bp4$ID %in% obsolete[,1],]
bp4 = bp4[order(bp4$p.adjust,decreasing=TRUE),]

#write.csv(bp4, "VSG_KEGG_trimed_down_prot_a.csv")

rk = character()
for(i in 1:(nrow(bp4)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp4)){
    pk = c(pk,sum(unlist(strsplit(bp4[i,8], "/")) %in% unlist(strsplit(bp4[j,8], "/")))/length(unlist(strsplit(bp4[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp4 = bp4[rk1 == "ostavi",]
#write.csv(bp4, "VSG_KEGG_trimed_down_prot_b.csv")

###############################################################################
# TRF - up
###############################################################################
setwd('/Users/viktorian.miok/Documents/consultation/Katarina/Integration_Data/Proteomics_data/results/TRF')

sig.gene <- bitr(ProteomOutput$diffTRF[(ProteomOutput$diffTRF$log2FoldChange > 0) & (ProteomOutput$diffTRF$pvalue < 0.05), 1],
                 fromType="SYMBOL",
                 toType="ENTREZID",
                 OrgDb=org.Mm.eg.db
)
kk1 <- enrichKEGG(gene=sig.gene[,2], 
                  organism='mmu'
)
ego1 <- enrichGO(gene=sig.gene[,2], 
                 OrgDb=org.Mm.eg.db, 
                 ont="BP", 
                 pAdjustMethod="BH",
                 readable=TRUE
)
ego2 <- enrichGO(gene=sig.gene[,2], 
                 OrgDb=org.Mm.eg.db,
                 ont="CC", 
                 pAdjustMethod="BH", 
                 readable=TRUE
)
ego3 <- enrichGO(gene=sig.gene[,2], 
                 OrgDb=org.Mm.eg.db, 
                 ont="MF", 
                 pAdjustMethod="BH",
                 readable=TRUE
)

bp1 = ego1@result
#write.csv(bp1, "TRF_GO_BP_up_prot.csv")

bp1 = bp1[bp1$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp1[,4]))
bp1 = bp1[(lk > 15) & (lk < 500),]
bp1 = bp1[bp1$Count > 4,]
bp1 = bp1[!bp1$ID %in% noliv[,1],]
bp1 = bp1[!bp1$ID %in% obsolete[,1],]
bp1 = bp1[order(bp1$p.adjust,decreasing=TRUE),]
#write.csv(bp1, "TRF_GO_BP_trimed_up_prot_a.csv")

rk = character()
for(i in 1:(nrow(bp1)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp1)){
    pk = c(pk,sum(unlist(strsplit(bp1[i,8], "/")) %in% unlist(strsplit(bp1[j,8], "/")))/length(unlist(strsplit(bp1[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp1 = bp1[rk1 == "ostavi",]
#write.csv(bp1, "TRF_GO_BP_trimed_up_prot_b.csv")

bp2 = ego2@result
#write.csv(bp2, "TRF_GO_CC_up_prot.csv")

bp2 = bp2[bp2$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp2[,4]))
bp2 = bp2[(lk > 15) & (lk < 500),]
bp2 = bp2[bp2$Count > 4,]
bp2 = bp2[!bp2$ID %in% noliv[,1],]
bp2 = bp2[!bp2$ID %in% obsolete[,1],]
bp2 = bp2[order(bp2$p.adjust,decreasing=TRUE),]
#write.csv(bp2, "TRF_GO_CC_trimed_up_prot_a.csv")

rk = character()
for(i in 1:(nrow(bp2)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp2)){
    pk = c(pk,sum(unlist(strsplit(bp2[i,8], "/")) %in% unlist(strsplit(bp2[j,8], "/")))/length(unlist(strsplit(bp2[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp2 = bp2[rk1 == "ostavi",]
#write.csv(bp2, "TRF_GO_CC_trimed_up_prot_b.csv")

bp3 = ego3@result
#write.csv(bp3, "TRF_GO_MF_up_prot.csv")

bp3 = bp3[bp3$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp3[,4]))
bp3 = bp3[(lk > 15) & (lk < 500),]
bp3 = bp3[bp3$Count > 4,]
bp3 = bp3[!bp3$ID %in% noliv[,1],]
bp3 = bp3[!bp3$ID %in% obsolete[,1],]
bp3 = bp3[order(bp3$p.adjust,decreasing=TRUE),]
#write.csv(bp3, "TRF_GO_MF_trimed_up_prot_a.csv")

rk = character()
for(i in 1:(nrow(bp3)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp3)){
    pk = c(pk,sum(unlist(strsplit(bp3[i,8], "/")) %in% unlist(strsplit(bp3[j,8], "/")))/length(unlist(strsplit(bp3[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp3 = bp3[rk1 == "ostavi",]
#write.csv(bp3, "TRF_GO_MF_trimed_up_prot_b.csv")

bp4 = kk1@result
for(i in 1:nrow(bp4)){
  bp4[i,8] <- paste(suppressMessages(bitr(unlist(strsplit(bp4[i,8], "/")),
                                          fromType="ENTREZID",
                                          toType="SYMBOL",
                                          OrgDb=org.Mm.eg.db)[,2]),
                    collapse="/"
  )
}

#write.csv(bp4, "TRF_KEGG_up_prot.csv")

bp4 = bp4[bp4$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp4[,4]))
bp4 = bp4[(lk > 15) & (lk < 500),]
bp4 = bp4[bp4$Count > 4,]
bp4 = bp4[!bp4$ID %in% noliv[,1],]
bp4 = bp4[!bp4$ID %in% obsolete[,1],]
bp4 = bp4[order(bp4$p.adjust,decreasing=TRUE),]

#write.csv(bp4, "TRF_KEGG_trimed_up_prot_a.csv")

rk = character()
for(i in 1:(nrow(bp4)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp4)){
    pk = c(pk,sum(unlist(strsplit(bp4[i,8], "/")) %in% unlist(strsplit(bp4[j,8], "/")))/length(unlist(strsplit(bp4[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp4 = bp4[rk1 == "ostavi",]
#write.csv(bp4, "TRF_KEGG_trimed_up_prot_b.csv")

###############################################################################
# TRF - down
###############################################################################
sig.gene <- bitr(ProteomOutput$diffTRF[(ProteomOutput$diffTRF$log2FoldChange < 0) & (ProteomOutput$diffTRF$pvalue < 0.05), 1],
                 fromType="SYMBOL",
                 toType="ENTREZID",
                 OrgDb=org.Mm.eg.db
)
kk1 <- enrichKEGG(gene=sig.gene[,2],
                  organism='mmu'
)
ego1 <- enrichGO(gene=sig.gene[,2], 
                 OrgDb=org.Mm.eg.db, 
                 ont="BP",
                 pAdjustMethod="BH",
                 readable=TRUE
)
ego2 <- enrichGO(gene=sig.gene[,2],
                 OrgDb=org.Mm.eg.db,
                 ont="CC",
                 pAdjustMethod="BH",
                 readable=TRUE
)
ego3 <- enrichGO(gene=sig.gene[,2],
                 OrgDb=org.Mm.eg.db, 
                 ont="MF", 
                 pAdjustMethod="BH", 
                 readable=TRUE
)
bp1 = ego1@result
#write.csv(bp1, "TRF_GO_BP_down_prot.csv")

bp1 = bp1[bp1$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp1[,4]))
bp1 = bp1[(lk > 15) & (lk < 500),]
bp1 = bp1[bp1$Count > 4,]
bp1 = bp1[!bp1$ID %in% noliv[,1],]
bp1 = bp1[!bp1$ID %in% obsolete[,1],]
bp1 = bp1[order(bp1$p.adjust,decreasing=TRUE),]# (as.numeric(sub("\\/.*", "", bp1[,4]))
#write.csv(bp1, "TRF_GO_BP_trimed_down_prot_a.csv")

rk = character()
for(i in 1:(nrow(bp1)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp1)){
    pk = c(pk,sum(unlist(strsplit(bp1[i,8], "/")) %in% unlist(strsplit(bp1[j,8], "/")))/length(unlist(strsplit(bp1[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp1 = bp1[rk1 == "ostavi",]
#write.csv(bp1, "TRF_GO_BP_trimed_down_prot_b.csv")

bp2 = ego2@result
#write.csv(bp2, "TRF_GO_CC_down_prot.csv")

bp2 = bp2[bp2$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp2[,4]))
bp2 = bp2[(lk > 15) & (lk < 500),]
bp2 = bp2[bp2$Count > 4,]
bp2 = bp2[!bp2$ID %in% noliv[,1],]
bp2 = bp2[!bp2$ID %in% obsolete[,1],]
bp2 = bp2[order(bp2$p.adjust,decreasing=TRUE),]
#write.csv(bp2, "TRF_GO_CC_trimed_down_prot_a.csv")

rk = character()
for(i in 1:(nrow(bp2)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp2)){
    pk = c(pk,sum(unlist(strsplit(bp2[i,8], "/")) %in% unlist(strsplit(bp2[j,8], "/")))/length(unlist(strsplit(bp2[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp2 = bp2[rk1 == "ostavi",]
#write.csv(bp2, "TRF_GO_CC_trimed_down_prot_b.csv")

bp3 = ego3@result
#write.csv(bp3, "TRF_GO_MF_down_prot.csv")

bp3 = bp3[bp3$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp3[,4]))
bp3 = bp3[(lk > 15) & (lk < 500),]
bp3 = bp3[bp3$Count > 4,]
bp3 = bp3[!bp3$ID %in% noliv[,1],]
bp3 = bp3[!bp3$ID %in% obsolete[,1],]
bp3 = bp3[order(bp3$p.adjust,decreasing=TRUE),]
#write.csv(bp3, "TRF_GO_MF_trimed_down_prot_a.csv")

rk = character()
for(i in 1:(nrow(bp3)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp3)){
    pk = c(pk,sum(unlist(strsplit(bp3[i,8], "/")) %in% unlist(strsplit(bp3[j,8], "/")))/length(unlist(strsplit(bp3[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp3 = bp3[rk1 == "ostavi",]
#write.csv(bp3, "TRF_GO_MF_trimed_down_prot_b.csv")

bp4 = kk1@result
for(i in 1:nrow(bp4)){
  bp4[i,8] <- paste(suppressMessages(bitr(unlist(strsplit(bp4[i,8], "/")),
                                          fromType="ENTREZID",
                                          toType="SYMBOL",
                                          OrgDb=org.Mm.eg.db)[,2]),
                    collapse="/"
)
}
#write.csv(bp4, "TRF_KEGG_down_prot.csv")

bp4 = bp4[bp4$p.adjust < 0.05,]
lk = as.numeric(sub("\\/.*", "", bp4[,4]))
bp4 = bp4[(lk > 15) & (lk < 500),]
bp4 = bp4[bp4$Count > 4,]
bp4 = bp4[!bp4$ID %in% noliv[,1],]
bp4 = bp4[!bp4$ID %in% obsolete[,1],]
bp4 = bp4[order(bp4$p.adjust,decreasing = TRUE),]

#write.csv(bp4, "TRF_KEGG_trimed_down_prot_a.csv")

rk = character()
for(i in 1:(nrow(bp4)-1)){
  pk = numeric()
  for(j in (i+1):nrow(bp4)){
    pk = c(pk,sum(unlist(strsplit(bp4[i,8], "/")) %in% unlist(strsplit(bp4[j,8], "/")))/length(unlist(strsplit(bp4[i,8], "/"))))
  }
  if (sum(pk > 0.4) > 0) {
    rk[i] = "brisi"
  } else rk[i] = "ostavi"
}
rk1 = c(rk, "ostavi")
bp4 = bp4[rk1 == "ostavi",]
#write.csv(bp4, "TRF_KEGG_trimed_down_prot_b.csv")




