rm(list=ls())

# load libraries
library("gridExtra")
pkgs=c("gridExtra",
        "DESeq2",
        "ggplot2",
        "gplots",
        "RColorBrewer",
        "clusterProfiler",
        "org.Mm.eg.db",
        "stringdist",
        "dplyr") # package names
inst=suppressMessages(lapply(pkgs, 
                             library,
                             character.only=TRUE)) # load them

setwd("~/Documents/consultation/Katarina/Integration_Data/Integration/data/")
load("Initial_Analysis.RData")

sampleName <- read.table("Sample_names.txt",
                         header=T, 
                         sep='\t'
)
colnames(expression.matrix) = as.character(sampleName[,2])
options(repr.plot.width=10, repr.plot.height=10)

# Obsolete terms
noliv <- read.csv("200828 non-liver related terms.csv",
                  header=FALSE
)
obsolete <- read.csv("200728_obsolete_GO_terms.csv",
                     header=FALSE
)
# Load functions
setwd("~/Documents/consultation/Katarina/Integration_Data/Integration/code/")
source("rna_prot_functions.R")
###############################################################################
#
#   TRF study
#
###############################################################################

trfdat = expression.matrix[,c(64:71,56:58,60:63)]

id = colnames(trfdat)
condition = as.factor(c(rep("TRFHFD",8), rep("ALHFD",7)))
metaData = data.frame(id, condition)

# differential expression analysis
TRF <- diff_expr_rna(data=trfdat,
                     mdata=metaData, 
                     design=~condition, 
                     filter=10
)
# volcano and PCA plot 
plot_fig_rna(deg_table=TRF$deg_tab,
             dds_object=TRF$dds, 
             col_volcano=c("grey", "dodgerblue3"),
             col_pca=c("dodgerblue4", "dodgerblue1"),
             title_volcano="TRF study - volcano plot",
             title_pca="TRF study - PCA plot",
             xlim=c(-5,5),
             ylim=c(0,20)
)
###############################################################################
#
#   VSG study
#
###############################################################################

vsgdat = expression.matrix[,c(75,77:81,83:89)]

id = colnames(vsgdat)
condition = as.factor(c(rep("VSGHFD",6), rep("PFHFD",7)))
metaData = data.frame(id, condition)

# differential expression analysis
VSG <- diff_expr_rna(data=vsgdat,
                     mdata=metaData, 
                     design=~condition, 
                     filter=10
)
# volcano and PCA plot 
plot_fig_rna(deg_table=VSG$deg_tab, 
             dds_object=VSG$dds,
             col_volcano=c("grey", "forestgreen"),
             col_pca=c("green4", "green1"), 
             title_volcano="TRF study - volcano plot",
             title_pca="TRF study - PCA plot",
             xlim=c(-5,5), 
             ylim=c(0,20)
)

###############################################################################
#
#   IF study
#
###############################################################################

ifdat = expression.matrix[,c(24,26:31,16:23)]  

id = colnames(ifdat)
condition = as.factor(c(rep("IFHFD",7), rep("ALHFD",8)))
metaData = data.frame(id, condition)

# differential expression analysis
IF <-  diff_expr_rna(data=ifdat,
                     mdata=metaData,
                     design=~ condition, 
                     filter=10
)

# volcano and PCA plot 
plot_fig_rna(deg_table=IF$deg_tab, 
             dds_object=IF$dds, col_volcano=c("grey", "firebrick"),
             col_pca=c("firebrick4", "firebrick1"),
             title_volcano="TRF study - volcano plot",
             title_pca="TRF study - PCA plot",
             xlim=c(-5,5), 
             ylim=c(0,20)
)


