

#remove all objects from the environment
rm(list=ls())

suppressMessages(library('readxl'))
suppressMessages(library('biomaRt'))

human <- useMart("ensembl",
                 dataset="hsapiens_gene_ensembl"
)
mouse <- useMart("ensembl", 
                 dataset="mmusculus_gene_ensembl"
)
# Set the directory
setwd('~/Documents/consultation/Katarina/mmu_to_hsa_geneID/')

###############################################################################
#
#       Up and down genes and proteins
#
###############################################################################

input = "211013 Up and down genes proteins.xlsx"

###############################
#   proteomics up-regulated
###############################
pr_up <- read_excel(input, 
                    sheet=1, 
                    skip=2, 
                    col_names=FALSE)

upIFprot <- getLDS(attributes=c("mgi_symbol"), 
                  filters="mgi_symbol", 
                  values=pr_up[[2]], 
                  mart=mouse,
                  attributesL=c("hgnc_symbol"),
                  martL=human, 
                  uniqueRows=TRUE
)
write.csv(upIFprot,
          "upIFprot.csv",
          row.names=FALSE
)

upTRFprot <- getLDS(attributes=c("mgi_symbol"), 
                   filters="mgi_symbol", 
                   values=pr_up[[4]], 
                   mart=mouse,
                   attributesL=c("hgnc_symbol"),
                   martL=human, 
                   uniqueRows=TRUE
)
write.csv(upTRFprot, "upTRFprot.csv", row.names=FALSE)

upVSGprot <- getLDS(attributes=c("mgi_symbol"), 
                   filters="mgi_symbol", 
                   values=pr_up[[6]], 
                   mart=mouse,
                   attributesL=c("hgnc_symbol"),
                   martL=human, 
                   uniqueRows=TRUE
)
write.csv(upVSGprot, 
          "upVSGprot.csv", 
          row.names=FALSE
)
###############################
#   proteomics down-regulated
###############################
pr_down <- read_excel(input,
                      sheet=2,
                      skip=2,
                      col_names=FALSE
)

downIFprot <- getLDS(attributes=c("mgi_symbol"), 
                    filters="mgi_symbol", 
                    values=pr_down[[2]], 
                    mart=mouse,
                    attributesL=c("hgnc_symbol"),
                    martL=human, 
                    uniqueRows=TRUE
)
write.csv(downIFprot, 
          "downIFprot.csv",
          row.names=FALSE
)

downTRFprot <- getLDS(attributes=c("mgi_symbol"), 
                     filters="mgi_symbol", 
                     values=pr_down[[4]], 
                     mart=mouse,
                     attributesL=c("hgnc_symbol"),
                     martL=human, 
                     uniqueRows=TRUE
)
write.csv(downTRFprot, 
          "downTRFprot.csv", 
          row.names=FALSE
)

downVSGprot <- getLDS(attributes=c("mgi_symbol"), 
                     filters="mgi_symbol", 
                     values=pr_down[[6]], 
                     mart=mouse,
                     attributesL=c("hgnc_symbol"),
                     martL=human, 
                     uniqueRows=TRUE
)
write.csv(downVSGprot, 
          "downVSGprot.csv",
          row.names=FALSE
)

########################################
#      transcriptomics  up-regulated
########################################
ge_up <- read_excel(input,
                    sheet=3, 
                    skip=2, 
                    col_names=FALSE)

upIFrna <- getLDS(attributes=c("mgi_symbol"), 
                 filters="mgi_symbol", 
                 values=ge_up[[1]], 
                 mart=mouse,
                 attributesL=c("hgnc_symbol"),
                 martL=human, 
                 uniqueRows=TRUE
)
write.csv(upIFrna, 
          "upIFrna.csv",
          row.names=FALSE
)

upTRFrna <- getLDS(attributes=c("mgi_symbol"), 
                  filters="mgi_symbol", 
                  values=ge_up[[3]], 
                  mart=mouse,
                  attributesL=c("hgnc_symbol"),
                  martL=human, 
                  uniqueRows=TRUE
)
write.csv(upTRFrna, 
          "upTRFrna.csv",
          row.names=FALSE
)

upVSGrna <- getLDS(attributes=c("mgi_symbol"), 
                  filters="mgi_symbol", 
                  values=ge_up[[5]], 
                  mart=mouse,
                  attributesL=c("hgnc_symbol"),
                  martL=human, 
                  uniqueRows=TRUE
)
write.csv(upVSGrna, 
          "upVSGrna.csv", 
          row.names=FALSE
)

########################################
#      transcriptomics  down-regulated
########################################
ge_down <- read_excel(input, sheet=4, skip=2, col_names=F)

downIFrna <- getLDS(attributes=c("mgi_symbol"), 
                   filters="mgi_symbol", 
                   values=ge_down[[1]], 
                   mart=mouse,
                   attributesL=c("hgnc_symbol"),
                   martL=human, 
                   uniqueRows=TRUE
)
write.csv(downIFrna, 
          "downIFrna.csv",
          row.names=FALSE
)

downTRFrna <- getLDS(attributes=c("mgi_symbol"), 
                    filters="mgi_symbol", 
                    values=ge_down[[3]], 
                    mart=mouse,
                    attributesL=c("hgnc_symbol"),
                    martL=human, 
                    uniqueRows=TRUE
)
write.csv(downTRFrna, 
          "downTRFrna.csv", 
          row.names=FALSE
)

downVSGrna <- getLDS(attributes=c("mgi_symbol"), 
                    filters="mgi_symbol", 
                    values=ge_down[[5]], 
                    mart=mouse,
                    attributesL=c("hgnc_symbol"),
                    martL=human, 
                    uniqueRows=TRUE
)
write.csv(downVSGrna,
          "downVSGrna.csv",
          row.names=FALSE
)

###############################################################################
#
#       All genes and proteins
#
###############################################################################

input = "Full gene_protein lists.xlsx"

###########################
#      proteomics
###########################
pr <- read_excel(input,
                 sheet=1, 
                 skip=3,
                 col_names=FALSE
)

aIFprot <- getLDS(attributes=c("mgi_symbol"), 
                 filters="mgi_symbol", 
                 values=pr[[2]], 
                 mart=mouse,
                 attributesL=c("hgnc_symbol"),
                 martL=human, 
                 uniqueRows=TRUE
)
write.csv(aIFprot, 
          "aIFprot.csv",
          row.names=FALSE
)

aTRFprot <- getLDS(attributes=c("mgi_symbol"), 
                  filters="mgi_symbol", 
                  values=pr[[4]], 
                  mart=mouse,
                  attributesL=c("hgnc_symbol"),
                  martL=human, 
                  uniqueRows=TRUE
)
write.csv(aTRFprot, 
          "aTRFprot.csv", 
          row.names=FALSE
)

aVSGprot <- getLDS(attributes=c("mgi_symbol"), 
                  filters="mgi_symbol", 
                  values=pr[[6]], 
                  mart=mouse,
                  attributesL=c("hgnc_symbol"),
                  martL=human, 
                  uniqueRows=TRUE
)
write.csv(aVSGprot, 
          "aVSGprot.csv",
          row.names=FALSE
)

dim(aVSGrna)
head(aVSGrna)

###########################
#      transcriptomics
###########################
ge <- read_excel(input,
                 sheet=2, 
                 skip=3, 
                 col_names=FALSE
)

aIFrna <- getLDS(attributes=c("mgi_symbol"), 
                 filters="mgi_symbol", 
                 values=ge[[1]], 
                 mart=mouse,
                 attributesL=c("hgnc_symbol"),
                 martL=human, 
                 uniqueRows=TRUE
)
write.csv(aIFrna,
          "aIFrna.csv", 
          row.names=FALSE
)

aTRFrna <- getLDS(attributes=c("mgi_symbol"), 
                  filters="mgi_symbol", 
                  values=ge[[3]], 
                  mart=mouse,
                  attributesL=c("hgnc_symbol"),
                  martL=human, 
                  uniqueRows=TRUE
)
write.csv(aTRFrna,
          "aTRFrna.csv", 
          row.names=FALSE
)

aVSGrna <- getLDS(attributes=c("mgi_symbol"), 
                  filters="mgi_symbol", 
                  values=ge[[5]], 
                  mart=mouse,
                  attributesL=c("hgnc_symbol"),
                  martL=human, 
                  uniqueRows=TRUE
)
write.csv(aVSGrna, 
          "aVSGrna.csv",
          row.names=FALSE
)

###############################################################################
#
#       Differential gene expressioin
#
###############################################################################

input = "210314 Proteomics_rnaseq differentially expressed.xlsx"

dgeIFprot <- getLDS(attributes=c("mgi_symbol"), 
                    filters="mgi_symbol", 
                    values=as.character(as.matrix(read_excel(input, sheet=1)[,2])),
                    mart=mouse,
                    attributesL=c("hgnc_symbol"),#,"chromosome_name", "start_position"),
                    martL=human, 
                    uniqueRows=TRUE
)
write.csv(dgeIFprot, 
          "dgeIFprot.csv",
          row.names=FALSE
)

dgeTRFprot <- getLDS(attributes=c("mgi_symbol"), 
                     filters="mgi_symbol", 
                     values=as.character(as.matrix(read_excel(input, sheet=2)[,2])),
                     mart=mouse,
                     attributesL=c("hgnc_symbol"),
                     martL=human, 
                     uniqueRows=TRUE
)
write.csv(dgeTRFprot, 
          "dgeTRFprot.csv", 
          row.names=FALSE
)

dgeVSGprot <- getLDS(attributes=c("mgi_symbol"), 
                     filters="mgi_symbol", 
                     values=as.character(as.matrix(read_excel(input, sheet=3)[,2])),
                     mart=mouse,
                     attributesL=c("hgnc_symbol"),
                     martL=human, 
                     uniqueRows=TRUE
)
write.csv(dgeVSGprot, 
          "dgeVSGprot.csv", 
          row.names=FALSE
)

dgeIFrna <- getLDS(attributes=c("mgi_symbol"), 
                   filters="mgi_symbol", 
                   values=as.character(as.matrix(read_excel(input, sheet=4)[,1])),
                   mart=mouse,
                   attributesL=c("hgnc_symbol"),
                   martL=human, 
                   uniqueRows=TRUE
)
write.csv(dgeIFrna, 
          "dgeIFrna.csv", 
          row.names=FALSE
)

dgeTRFrna <- getLDS(attributes=c("mgi_symbol"), 
                    filters="mgi_symbol", 
                    values=as.character(as.matrix(read_excel(input, sheet=5)[,1])),
                    mart=mouse,
                    attributesL=c("hgnc_symbol"),
                    martL=human, 
                    uniqueRows=TRUE
)
write.csv(dgeTRFrna, 
          "dgeTRFrna.csv", 
          row.names=FALSE
)

dgeVSGrna <- getLDS(attributes=c("mgi_symbol"), 
                    filters="mgi_symbol", 
                    values=as.character(as.matrix(read_excel(input, sheet=6)[,1])),
                    mart=mouse,
                    attributesL=c("hgnc_symbol"),
                    martL=human, 
                    uniqueRows=TRUE
)
write.csv(dgeVSGrna, 
          "dgeVSGrna.csv", 
          row.names=FALSE
)

###############################################################################
#
#       Gene enrichment analysis
#
###############################################################################

input1 = "210314 Proteomics_rnaseq from enrichment analysis.xlsx"

geaIFprot <- getLDS(attributes=c("mgi_symbol"), 
                    filters="mgi_symbol", 
                    values=as.character(as.matrix(read_excel(input1, sheet=1)[,1])),
                    mart=mouse,
                    attributesL=c("hgnc_symbol"),#,"chromosome_name", "start_position"),
                    martL=human, 
                    uniqueRows=TRUE
)
write.csv(geaIFprot,
          "geaIFprot.csv",
          row.names=FALSE
)

geaTRFprot <- getLDS(attributes=c("mgi_symbol"), 
                     filters="mgi_symbol", 
                     values=as.character(as.matrix(read_excel(input1, sheet=2)[,1])),
                     mart=mouse,
                     attributesL=c("hgnc_symbol"),
                     martL=human, 
                     uniqueRows=TRUE
)
write.csv(geaTRFprot, 
          "geaTRFprot.csv", 
          row.names=FALSE)

geaVSGprot <- getLDS(attributes=c("mgi_symbol"), 
                     filters="mgi_symbol", 
                     values=as.character(as.matrix(read_excel(input1, sheet=3)[,1])),
                     mart=mouse,
                     attributesL=c("hgnc_symbol"),
                     martL=human, 
                     uniqueRows=TRUE
)
write.csv(geaVSGprot, 
          "geaVSGprot.csv",
          row.names=FALSE
)

geaIFrna <- getLDS(attributes=c("mgi_symbol"), 
                   filters="mgi_symbol", 
                   values=as.character(as.matrix(read_excel(input1, sheet=4)[,1])),
                   mart=mouse,
                   attributesL=c("hgnc_symbol"),
                   martL=human, 
                   uniqueRows=TRUE
)
write.csv(geaIFrna,
          "geaIFrna.csv",
          row.names=FALSE
)

geaTRFrna <- getLDS(attributes=c("mgi_symbol"), 
                    filters="mgi_symbol", 
                    values=as.character(as.matrix(read_excel(input1, sheet=5)[,1])),
                    mart=mouse,
                    attributesL=c("hgnc_symbol"),
                    martL=human, 
                    uniqueRows=TRUE
)
write.csv(geaTRFrna, 
          "geaTRFrna.csv", 
          row.names=FALSE
)

geaVSGrna <- getLDS(attributes=c("mgi_symbol"), 
                    filters="mgi_symbol", 
                    values=as.character(as.matrix(read_excel(input1, sheet=6)[,1])),
                    mart=mouse,
                    attributesL=c("hgnc_symbol"),
                    martL=human, 
                    uniqueRows=TRUE
)
write.csv(geaVSGrna, 
          "geaVSGrna.csv",
          row.names=FALSE
)




