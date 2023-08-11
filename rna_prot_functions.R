
###############################################################################
#
#       Function for differential gene expression analysis of transcriptomics data
#   
###############################################################################

diff_expr_rna <- function(data, mdata, design, filter){
  dds <- DESeqDataSetFromMatrix(countData=data,
                                colData=mdata, 
                                design=design
  )
  featureData = data.frame(gene=rownames(dds))
  mcols(dds) = DataFrame(mcols(dds), 
                         featureData
  )

  design(dds)=design 
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds, 
                             fit = 'parametric'
  )
  
  keep = rowSums(counts(dds)) >= filter
  dds = dds[keep,]
  
  dds <- DESeq(dds)
  res_deg = results(dds)
  
  res = res_deg[order(res_deg$pvalue),]
  deg = as.data.frame(dplyr::mutate(as.data.frame(res), 
                                    significant=ifelse(res$pvalue < 0.05, "pvalue<0.05", "not signif.")), 
                      row.names=rownames(res))
  deg = deg[(!is.na(deg$pvalue)),]
  return(list(deg = deg, 
              dds = dds)) 
}

###############################################################################
#
#       Function for generating volcano and PCA plot of transcriptomics data
#   
###############################################################################

plot_fig_rna <- function(deg_table, dds_object, col_volcano, col_pca, xlim, ylim){
  p1 <- ggplot(deg_table, aes(log2FoldChange, -log10(pvalue))) +
               geom_point(aes(col=significant)) +
               scale_color_manual(values=col_volcano) +
               guides(col=guide_legend(nrow=2)) + 
               xlim(xlim) + 
               ylim(ylim)
  
  vsd <- assay(vst(dds_object, 
                   fitType="local",
                   blind = TRUE)
  )
  df_pca <- prcomp(t(vsd))
  df_out = as.data.frame(df_pca$x)
  
  p2 <- ggplot(df_out,aes(x=PC1,y=PC2,color=condition)) + 
               scale_color_manual(values = col_pca) +
               geom_point(size = 5)
  
  grid.arrange(p1, p2, ncol=2)
}
