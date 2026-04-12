## Load packages

suppressMessages(library(DESeq2))
library(ggplot2)
library(tidyverse)
library(AnnotationDbi)
library(EnsDb.Hsapiens.v86)
library(biomaRt)
############---------------------------------------------------------------------



## HELPER FXN

ggsave_png.svg = function(plot.obj, h=5, w =5, folder, name){
  
  ## save png
  ggsave(plot = plot.obj, filename = paste0(folder, "/", name, ".png"),
         height = h, width = w)
  ## save svg
  ggsave(plot = plot.obj, filename = paste0(folder, "/", name, ".svg"),
         height = h, width = w)
  
}


counts <- read.table("path/to/BryoJQ1_AllSamples.Rmatrix.txt",
                     sep = "\t",

coldata <- read.csv("path/to/BryoJQ1_coldata.csv",
                    row.names = 1, header = T, sep = ";")

## get rRNA genes to remove from counts

mart = useEnsembl('ensembl', dataset = "hsapiens_gene_ensembl")
rrna_genes = biomaRt::getBM(values="rRNA", 
                            filters="biotype", 
                            attributes=c("ensembl_gene_id", "external_gene_name", "gene_biotype"), 
                            mart = mart)


## filter out all counts aligning to rrna genes
counts_noRRNA= counts[!(rownames(counts) %in% rrna_genes$ensembl_gene_id),]

## create grouping factor 

coldata$group = paste0(coldata$clone, "_", coldata$treatment)
coldata$replicate = paste0("R",coldata$replicate)
coldata

## create DESEQ2 dds #####

## exchange colnames for rownames of coldata
colnames(counts_noRRNA) = rownames(coldata)
counts_noRRNA

## rearrange 
counts_noRRNA2 = cbind(counts_noRRNA[,1:20], counts_noRRNA[22], counts_noRRNA[21], counts_noRRNA[,23:24])
counts_noRRNA2

coldata2 = rbind(coldata[1:20,], coldata[22,], coldata[21,], coldata[23:24,])
coldata2

# by clone_treat group ####
dds <- DESeqDataSetFromMatrix(countData = counts_noRRNA2,
                              colData = coldata2,
                              design = ~ 0 + group + replicate)
dds<- DESeq(dds)
## remove genes with fewer than 10 summed counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
counts(dds)
## save dds

saveRDS(dds, file= paste0("path/to/",
                          "dds_clone_treat_group.rds"))



## plot PCA
vst_dds <- vst(dds)

data <- plotPCA(vst_dds, intgroup = c( "group"), returnData=TRUE)
data$clone = gsub("(.*)_(.*)", replacement = "\\1", x = data$group)
data$treatment = gsub("(.*)_(.*)", replacement = "\\2", x = data$group)
data

percentVar <- round(100 * attr(data, "percentVar"))

pca <- ggplot(data, aes(PC1, PC2,  shape = clone, colour = treatment)) + 
  geom_point(size=3) +
  scale_colour_manual(values = c("mock"= "grey",
                                 "bryo" = "#E5625E",
                                 "jq1" = "#0C98E9",
                                 "bryo.jq1" = "#9B4CAC"))+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_bw()

pca














######################################################################################################

# by treatonly group ####
dds2 <- DESeqDataSetFromMatrix(countData = counts_noRRNA,
                               colData = coldata,
                               design = ~ 0 + treatment + replicate)
dds2<- DESeq(dds2)
## remove genes with fewer than 10 summed counts
keep <- rowSums(counts(dds2)) >= 10
dds2 <- dds2[keep,]

## save dds

saveRDS(dds2, file= paste0("C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/HIV_LAT.CLONE/LRA_testing/JLEG2_Bryo_JQ1_synergy/RNAseq_ATACseq COMBO/rnaseq",
                           "dds_treatonly_group.rds"))



## plot PCA
vst_dds2 <- vst(dds2)

data <- plotPCA(vst_dds2, intgroup = c( "group"), returnData=TRUE)
#data$clone = gsub("(.*)_(.*)", replacement = "\\1", x = data$group)
data$treatment = gsub("(.*)_(.*)", replacement = "\\2", x = data$group)
data

percentVar <- round(100 * attr(data, "percentVar"))

pca <- ggplot(data, aes(PC1, PC2, colour = treatment)) + 
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_classic()

pca






















####################################

data$treat <- gsub("(.*)_.*", "\\1", data$group)
data$infection <- gsub(".*_(.*)", "\\1", data$group)

data

res <- results(dds, contrast=c(1,0,0,0,-1,0,0,0,0,0))
res2 = as.data.frame(res) %>%
  filter(padj < 0.05)
res2


res <- results(dds, contrast=c(0,0,0,1,0,0,0,-1,0,0))
res2 = as.data.frame(res) %>%
  filter(padj < 0.01)
res2


res <- results(dds, contrast=c(1,0,0,-1,1,0,0,-1,0,0))
res2 = as.data.frame(res) %>%
  filter(padj < 0.05)


## 2.1 bryo v mock
res <- results(dds, contrast=c(1,0,0,-1,0,0,0,0,0,0))
res2 = as.data.frame(res) %>%
  dplyr::filter(padj < 0.01)

res2$ensID = rownames(res2)

res2$gene <- mapIds(EnsDb.Hsapiens.v86, keys= res2$ensID, column="SYMBOL", keytype="GENEID", multiVals="first")
head(res2$gene, n=100)

res3 = as.data.frame(res2)

## 2.un1 bryo v mock
res <- results(dds, contrast=c(0,0,0,0,1,0,0,-1,0,0))
res2 = as.data.frame(res) %>%
  dplyr::filter(padj < 0.01)

res2$ensID = rownames(res2)

res2$gene <- mapIds(EnsDb.Hsapiens.v86, keys= res2$ensID, column="SYMBOL", keytype="GENEID", multiVals="first")
head(res2$gene, n=100)

res4 = as.data.frame(res2)


test = setdiff(x = res3$gene, res4$gene)


coldata
res2$gene[res2$gene == "SERINC5"]
res2


## 2.1 jq1 v mock
res <- results(dds, contrast=c(0,0,1,-1,0,0,0,0,0,0))
res <- results(dds, contrast=c(0,0,0,0,0,0,1,-1,0,0))
res2 = as.data.frame(res) %>%
  dplyr::filter(padj < 0.01)

res2$ensID = rownames(res2)

res2$gene <- mapIds(EnsDb.Hsapiens.v86, keys= res2$ensID, column="SYMBOL", keytype="GENEID", multiVals="first")
head(res2$gene, n=100)

#res6 = as.data.frame(res2)
res7 = as.data.frame(res2)
res5 = res3[res3$gene %in% test,]
res5
test

## 2.1 bryojq1 v bryo
res <- results(dds, contrast=c(-1,1,0,0,0,0,0,0,0,0))
res2 = as.data.frame(res) %>%
  dplyr::filter(padj < 0.05)

res2$ensID = rownames(res2)

res2$gene <- mapIds(EnsDb.Hsapiens.v86, keys= res2$ensID, column="SYMBOL", keytype="GENEID", multiVals="first")
head(res2$gene, n=100)

res8 = as.data.frame(res2)
