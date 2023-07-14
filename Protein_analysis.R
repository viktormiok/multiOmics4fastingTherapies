rm(list=ls())

# load libraries
pkgs = c("DEP",
         "GenomicRanges",
         "SummarizedExperiment",
         "gplots",
         "dplyr",
         "EnhancedVolcano",
         "gdata",
         "clusterProfiler",
         "org.Mm.eg.db",
         "VennDiagram",
         "RColorBrewer",
         "readxl") # package names
inst <- suppressMessages(lapply(pkgs,
                                library, 
                                character.only=TRUE)
) # load them

###############################################################################
#
#     IF Study
#
###############################################################################
setwd('~/Documents/consultation/Katarina/Integration_Data/Proteomics_data')
dat <- read_excel('IF_study.xlsx',
                  sheet=3
)
colors = colorRampPalette(rev(brewer.pal(9, "Spectral")))(255)
dat = dat[1:2215,]

data_unique <- make_unique(dat,
                           "PG.Genes",
                           "PG.ProteinDescriptions",
                           delim="/t"
)
a_columns = c(grep(c("IFHFD_"), colnames(data_unique)),grep(c("ALHFD_"), colnames(data_unique)))
# removed sample
a_columns = a_columns[-4]

label = colnames(data_unique)[a_columns]
condition = c("IFHFD","IFHFD","IFHFD","IFHFD","IFHFD","IFHFD","IFHFD",
              "ALHFD","ALHFD","ALHFD","ALHFD","ALHFD","ALHFD","ALHFD")

replicate = c(1:7,1:7)
experimental_design = data.frame(label, condition, replicate)  
experimental_design[,1] = as.character(experimental_design[,1])
experimental_design[,2] = as.character(experimental_design[,2])
experimental_design[,3] = as.character(experimental_design[,3])

data_se <- make_se(data_unique,
                   a_columns,
                   experimental_design
)
assay(data_se)[is.nan(assay(data_se))] = NA

# Plot a barplot of the protein identification overlap between samples
plot_frequency(data_se)

# Filter for proteins that are identified in 2 out of 3 replicates of at least one condition
data_filt <- filter_missval(data_se,
                            thr=1
)
# Plot a barplot of the number of identified proteins per samples
plot_numbers(data_filt)

# Plot a barplot of the protein identification overlap between samples
plot_coverage(data_filt)

# Normalize the data
data_norm <- normalize_vsn(data_filt)

# Visualize normalization by boxplots for all samples before and after normalization
plot_normalization(data_filt,
                   data_norm
)

# Plot a heatmap of proteins with missing values
plot_missval(data_filt)

# Plot intensity distributions and cumulative fraction of proteins with and without missing values
plot_detect(data_filt)

# Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
data_imp <- impute(data_norm,
                   fun="MinProb",
                   q=0.01
)
# Plot intensity distributions before and after imputation
plot_imputation(data_norm,
                data_imp
)

data_diff_manual <- test_diff(data_imp,
                              type="manual", 
                              test=c("IFHFD_vs_ALHFD")
)

# Denote significant proteins based on user defined cutoffs
depIF <- add_rejections(data_diff_manual, 
                        alpha=0.05, 
                        lfc=log2(2)
)
res_ALHFD_vs_IFHFD <- get_results(depIF)

# Add entrezID column
#sig.gene <- bitr(res_ALHFD_vs_IFHFD[,1],
#                 fromType="SYMBOL",
#                 toType="ENTREZID",
#                 OrgDb='org.Mm.eg.db'
#)
#colnames(sig.gene) = c("name","ENTREZID")
#res_ALHFD_vs_IFHFD_entrez <- merge(sig.gene,
#                                   res_ALHFD_vs_IFHFD,
#                                   by="name"
#)

#write.table(res_ALHFD_vs_IFHFD, 
#            file="IFstudy_IFHFDvsALHFD.csv",
#            sep=",",
#            col.names=NA
#)
#write.table(res_ALHFD_vs_IFHFD_entrez,
#            file="IFstudy_IFHFDvsALHFD_entrez.csv",
#            sep=",",
#            col.names=NA
#)

colnames(res_ALHFD_vs_IFHFD) = c("name","ID","pvalue","padj","sig","sig1",
                                 "log2FoldChange","AL_centered","IF_centered")

res = res_ALHFD_vs_IFHFD[order(res_ALHFD_vs_IFHFD$pvalue),]
resultsIF = as.data.frame(dplyr::mutate(as.data.frame(res),
                                        significant=ifelse(res$pvalue < 0.05, "pvalue<0.05", "not signif.")),
                        row.names=rownames(res)
)

p <- ggplot2::ggplot(resultsIF, ggplot2::aes(log2FoldChange, -log10(pvalue))) +
     ggplot2::geom_point(ggplot2::aes(col=significant)) +
     ggplot2::scale_color_manual(values=c("grey", "firebrick")) +
     guides(col=guide_legend(nrow=2)) +
     ggplot2::ggtitle("") +
     xlim(-2.5,2.5) +
     ylim(0,7)
p

df_pca = prcomp(t(assay(data_imp)))
df_out = as.data.frame(df_pca$x)

p <- ggplot(df_out,aes(x=PC1,y=PC2,color=condition)) + 
     ggplot2::scale_color_manual(values = c("firebrick4", "firebrick1")) +
     geom_point(size = 5)
p

###############################################################################
#
#     TRF Study
#
###############################################################################

dat <- read_excel('TRF_study.xlsx',
                  sheet=3
)
colors = colorRampPalette(rev(brewer.pal(9, "Spectral")))(255)

data_unique <- make_unique(dat, 
                           "PG.Genes",
                           "PG.ProteinDescriptions", 
                           delim="/t"
)
a_columns = c(grep(c("TRFHFD_"), colnames(data_unique)),
              grep(c("ALHFD_"), colnames(data_unique)))
# removed sample
a_columns = a_columns[-10]

label = colnames(data_unique)[a_columns]
condition = c("TRFHFD","TRFHFD","TRFHFD","TRFHFD","TRFHFD","TRFHFD",
              "ALHFD","ALHFD","ALHFD","ALHFD","ALHFD")

replicate = c(1:6,1:5)
experimental_design = data.frame(label, condition, replicate)  
experimental_design[,1] = as.character(experimental_design[,1])
experimental_design[,2] = as.character(experimental_design[,2])
experimental_design[,3] = as.character(experimental_design[,3])

data_se <- make_se(data_unique, 
                   a_columns,
                   experimental_design
)
assay(data_se)[is.nan(assay(data_se))] = NA

# Plot a barplot of the protein identification overlap between samples
plot_frequency(data_se)

# Filter for proteins that are identified in 2 out of 3 replicates of at least one condition
data_filt <- filter_missval(data_se,
                            thr = 1
)
# Plot a barplot of the number of identified proteins per samples
plot_numbers(data_filt)

# Plot a barplot of the protein identification overlap between samples
plot_coverage(data_filt)

# Normalize the data
data_norm <- normalize_vsn(data_filt)

# Visualize normalization by boxplots for all samples before and after normalization
plot_normalization(data_filt,
                   data_norm
)
# Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
data_imp <- impute(data_norm, 
                   fun="MinProb", 
                   q=0.01
)
# Plot intensity distributions before and after imputation
plot_imputation(data_norm, 
                data_imp
)

data_diff_manual <- test_diff(data_imp,
                              type="manual", 
                              test=c("TRFHFD_vs_ALHFD")
)
# Denote significant proteins based on user defined cutoffs
depTRF <- add_rejections(data_diff_manual,
                         alpha=0.05,
                         lfc=log2(2)
)
res_TRFHFD_vs_ALHFD = get_results(depTRF)

# Add entrezID column
#sig.gene <- bitr(res_TRFHFD_vs_ALHFD[,1], 
#                 fromType="SYMBOL",
#                 toType="ENTREZID",
#                 OrgDb='org.Mm.eg.db'
#)
#colnames(sig.gene) = c("name","ENTREZID")
#res_TRFHFD_vs_ALHFD_entrez <- merge(sig.gene,
#                                    res_TRFHFD_vs_ALHFD,
#                                    by="name"
#)

#write.table(res_TRFHFD_vs_ALHFD, 
#            file="TRFstudy_TRFHFDvsALHFD.csv", 
#            sep=",",
#            col.names=NA
#)
#write.table(res_TRFHFD_vs_ALHFD_entrez, 
#            file="TRFstudy_TRFHFDvsALHFD_entrez.csv",
#            sep=",",
#            col.names=NA
#)

colnames(res_TRFHFD_vs_ALHFD) = c("name","ID","pvalue","padj","sig","sig1",
                                  "log2FoldChange","AL_centered","IF_centered")

res = res_TRFHFD_vs_ALHFD[order(res_TRFHFD_vs_ALHFD$pvalue),]
resultsTRF = as.data.frame(dplyr::mutate(as.data.frame(res), 
                                         significant=ifelse(res$pvalue < 0.05, "pvalue<0.05", "not signif.")),
                            row.names=rownames(res))

p <- ggplot2::ggplot(resultsTRF, ggplot2::aes(log2FoldChange, -log10(pvalue))) +
     ggplot2::geom_point(ggplot2::aes(col=significant)) +
     ggplot2::scale_color_manual(values=c("grey", "dodgerblue3")) +
     guides(col = guide_legend(nrow=2)) + 
     ggplot2::ggtitle("") +
     xlim(-2.5,2.5) +
     ylim(0,7)
p

df_pca <- prcomp(t(assay(data_imp)))
df_out = as.data.frame(df_pca$x)

p <- ggplot(df_out, aes(x=PC1, y=PC2, color=condition)) +
     ggplot2::scale_color_manual(values = c("dodgerblue4", "dodgerblue1")) +
     geom_point(size = 5)
p

###############################################################################
#
#     VSG Study
#
###############################################################################

dat <- read_excel('191119 VSG study KLK13007_DirectDIA_nodecimal.xlsx',
                  sheet=2
)
data_unique <- make_unique(dat,
                           "PG.Genes",
                           "PG.ProteinDescriptions", 
                           delim="/t"
)
a_columns = c(grep(c("VSG HFD_"), colnames(data_unique)),
              grep(c("PF HFD_"), colnames(data_unique)))

# removed sample
a_columns = a_columns[-c(2,9)]

label = colnames(data_unique)[a_columns]
condition = c("VSGHFD","VSGHFD","VSGHFD","VSGHFD","VSGHFD","VSGHFD","VSGHFD",
              "PFHFD","PFHFD","PFHFD","PFHFD","PFHFD","PFHFD","PFHFD")

replicate = c(1:7,1:7)
experimental_design = data.frame(label, condition, replicate)  
experimental_design[,1] = as.character(experimental_design[,1])
experimental_design[,2] = as.character(experimental_design[,2])
experimental_design[,3] = as.character(experimental_design[,3])

#data_unique <- apply(data_unique[,c(7,9:14,42:48)], 2, as.integer)

data_unique[,7] = as.integer(data_unique[,7])
data_unique[,9] = as.integer(data_unique[,9])
data_unique[,10] = as.integer(data_unique[,10])
data_unique[,11] = as.integer(data_unique[,11])
data_unique[,12] = as.integer(data_unique[,12])
data_unique[,13] = as.integer(data_unique[,13])
data_unique[,14] = as.integer(data_unique[,14])

data_unique[,42] = as.integer(data_unique[,42])
data_unique[,43] = as.integer(data_unique[,43])
data_unique[,44] = as.integer(data_unique[,44])
data_unique[,45] = as.integer(data_unique[,45])
data_unique[,46] = as.integer(data_unique[,46])
data_unique[,47] = as.integer(data_unique[,47])
data_unique[,48] = as.integer(data_unique[,48])

data_se <- make_se(data_unique, 
                   a_columns, 
                   experimental_design
)

assay(data_se)[is.nan(assay(data_se))] = NA

# Plot a barplot of the protein identification overlap between samples
plot_frequency(data_se)

# Filter for proteins that are identified in 2 out of 3 replicates of at least one condition
data_filt <- filter_missval(data_se,
                            thr=1)

# Plot a barplot of the number of identified proteins per samples
plot_numbers(data_filt)

# Plot a barplot of the protein identification overlap between samples
plot_coverage(data_filt)

# Normalize the data
data_norm <- normalize_vsn(data_filt)

# Visualize normalization by boxplots for all samples before and after normalization
plot_normalization(data_filt, 
                   data_norm
)

# Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
data_imp <- impute(data_norm, 
                   fun="MinProb",
                   q=0.01
)
# Plot intensity distributions before and after imputation
plot_imputation(data_norm, 
                data_imp
)
data_diff_manual <- test_diff(data_imp,
                              type="manual", 
                              test=c("VSGHFD_vs_PFHFD")
)
# Denote significant proteins based on user defined cutoffs
depVGS <- add_rejections(data_diff_manual, 
                         alpha=0.05,
                         lfc=log2(2)
)
res_VSGHFD_vs_PFHFD = get_results(depVGS)

# Add entrezID column
sig.gene <- bitr(res_VSGHFD_vs_PFHFD[,1],
                 fromType="SYMBOL",
                 toType="ENTREZID",
                 OrgDb='org.Mm.eg.db'
)
colnames(sig.gene) = c("name","ENTREZID")
res_VSGHFD_vs_PFHFD_entrez <- merge(sig.gene,
                                    res_VSGHFD_vs_PFHFD,
                                    by="name"
)

write.table(res_VSGHFD_vs_PFHFD, 
            file="res_VSGHFD_vs_PFHFD.csv",
            sep=",",
            col.names=NA
)
write.table(res_VSGHFD_vs_PFHFD_entrez,
            file="res_VSGHFD_vs_PFHFD_entrez.csv",
            sep=",",
            col.names=NA
)

colnames(res_VSGHFD_vs_PFHFD) = c("name","ID","pvalue","padj","sig","sig1",
                                  "log2FoldChange","AL_centered","IF_centered")

res=res_VSGHFD_vs_PFHFD[order(res_VSGHFD_vs_PFHFD$pvalue),]
resultsVSG=as.data.frame(dplyr::mutate(as.data.frame(res),
                                       significant=ifelse(res$pvalue < 0.05, "pvalue<0.05", "not signif.")),
                         row.names=rownames(res))

p <- ggplot2::ggplot(resultsVSG, ggplot2::aes(log2FoldChange, -log10(pvalue))) +
     ggplot2::geom_point(ggplot2::aes(col=significant)) + 
     ggplot2::scale_color_manual(values=c("grey", "forestgreen")) +
     guides(col=guide_legend(nrow=2)) +
     ggplot2::ggtitle("") + 
     xlim(-2.5,2.5) + 
     ylim(0,7)
p

df_pca = prcomp(t(assay(data_imp)))
df_out = as.data.frame(df_pca$x)

p <- ggplot(df_out,aes(x=PC1,y=PC2,color=condition)) +
     ggplot2::scale_color_manual(values=c("green4", "green1")) + 
     geom_point(size=5)
p





