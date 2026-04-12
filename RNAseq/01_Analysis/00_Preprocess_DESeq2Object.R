## Load packages

suppressMessages(library(DESeq2))
library(ggplot2)
library(tidyverse)
library(AnnotationDbi)
library(EnsDb.Hsapiens.v86)
library(biomaRt)
#--------------------------------------------------------------------



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
                     row.names = 1, header = T)

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

