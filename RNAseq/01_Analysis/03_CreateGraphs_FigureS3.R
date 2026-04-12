## load packages

suppressMessages(library(DESeq2))
library(ggplot2)
library(tidyverse)
library(dplyr)
library(AnnotationDbi)
library(EnsDb.Hsapiens.v86)
library(biomaRt)
library(ggvenn)
library(pheatmap)
library(ggplotify)
library(patchwork)
library(heatmaply)
library(msigdbr)
library(ggupset)
library(ggprism)

#------------------------------------------------------------------------------

## read in differential gene expression results

deg_list = readRDS("res_list.rds")

## bind together all contrast results into one dataframe for graphing
res_mega = deg_list[c(1:3, 5:7)] %>%
  bind_rows(., .id = "contrast") 

## alter names of contrasts for better graphical representation

res_mega = res_mega %>%
  mutate(Clone = ifelse(grepl("jleg2.1", contrast), "JLEG2.1", "Uninfect")) %>%
  mutate(Treat = ifelse(grepl("bryo_v_mock", contrast), "Bryo",
                        ifelse(grepl("1.jq1_v_mock", contrast), "JQ1", 
                               ifelse(grepl("bryo.jq1_v_mock", contrast), "BryoJQ1", "TEST"))))  %>%
  mutate(CloneTreat = paste0(Clone, ".", Treat))


#### Figure S3A ####

## Create heatmap of top 20 up and downregulated genes by contrast


### set palette
paletteLength <- 50
myColor <- colorRampPalette(c('#0000FF', "white", "darkorange"))(paletteLength)

myBreaks <- c(seq(-9, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(5/paletteLength, 9, length.out=floor(paletteLength/2)))

# top 20 upregulated genes
top20_up = res_mega %>%
  group_by(contrast) %>%
  arrange(-log2FoldChange, .by_group = T) %>%
  top_n(20, log2FoldChange) %>%
  pull(gene) %>%
  unique(.)
# top 20 downregulated genes
top20_down = res_mega %>%
  group_by(contrast) %>%
  arrange(log2FoldChange, .by_group = T) %>%
  top_n(-20, log2FoldChange) %>%
  pull(gene) %>%
  unique(.)

## make heatmaps for up and down

## up 
sub = res_mega %>%
  mutate(log2FoldChange = ifelse(padj > 0.05, 0, log2FoldChange)) %>%
  mutate(log2FoldChange = ifelse(is.na(log2FoldChange), 0, log2FoldChange)) %>%
  #dplyr::filter(!is.na(log2FoldChange)) %>%
  dplyr::filter( gene %in% top20_up)%>%
  select(gene, log2FoldChange, CloneTreat) %>%
  pivot_wider(names_from = CloneTreat, values_from = log2FoldChange, id_cols = gene, values_fn = function(x) mean(x)) %>%
  column_to_rownames(var = "gene") %>%
  relocate("Uninfect.Bryo", "Uninfect.JQ1", "Uninfect.BryoJQ1", 
           "JLEG2.1.Bryo", "JLEG2.1.JQ1", "JLEG2.1.BryoJQ1") %>%
  data.matrix()
sub


hmap.top20_up = pheatmap(sub, color=myColor, breaks=myBreaks,  cluster_cols = F, gaps_col = 3)#, 

## down
sub = res_mega %>%
  mutate(log2FoldChange = ifelse(padj > 0.05, 0, log2FoldChange)) %>%
  mutate(log2FoldChange = ifelse(is.na(log2FoldChange), 0, log2FoldChange)) %>%
  #dplyr::filter(!is.na(log2FoldChange)) %>%
  dplyr::filter( gene %in% top20_down)%>%
  select(gene, log2FoldChange, CloneTreat) %>%
  pivot_wider(names_from = CloneTreat, values_from = log2FoldChange, id_cols = gene, values_fn = function(x) mean(x)) %>%
  column_to_rownames(var = "gene") %>%
  relocate("Uninfect.Bryo", "Uninfect.JQ1", "Uninfect.BryoJQ1", 
           "JLEG2.1.Bryo", "JLEG2.1.JQ1", "JLEG2.1.BryoJQ1") %>%
  data.matrix()
sub


hmap.top20_down = pheatmap(sub, color=myColor, breaks=myBreaks,  cluster_cols = F, gaps_col = 3)


#### Figure S3B ####

## Create heatmap of IFN signalling gene panel

### set palette
paletteLength <- 50
myColor <- colorRampPalette(c('#0000FF', "white", "darkorange"))(paletteLength)

myBreaks <- c(seq(-5, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(5/paletteLength, 5, length.out=floor(paletteLength/2)))

## retrieve Reactome IFN signalling pathway
reacto_ifngsignal = msigdbr(species = "human", category = "C2", subcollection = "CP:REACTOME")%>%
  dplyr::filter(gs_exact_source == "R-HSA-877300") %>%
  pull(unique(gene_symbol))
reacto_ifngsignal

sub = res_mega %>%
  mutate(log2FoldChange = ifelse(padj > 0.05, 0, log2FoldChange)) %>%
  mutate(log2FoldChange = ifelse(is.na(log2FoldChange), 0, log2FoldChange)) %>%
  dplyr::filter( gene %in% reacto_ifngsignal)%>%
  select(gene, log2FoldChange, CloneTreat) %>%
  pivot_wider(names_from = CloneTreat, values_from = log2FoldChange, id_cols = gene, values_fn = function(x) mean(x)) %>%
  column_to_rownames(var = "gene") %>%
  relocate("Uninfect.Bryo", "Uninfect.JQ1", "Uninfect.BryoJQ1", 
           "JLEG2.1.Bryo", "JLEG2.1.JQ1", "JLEG2.1.BryoJQ1") %>%
  data.matrix()

hmap.ifngsignal = pheatmap(sub, color=myColor, breaks=myBreaks,  cluster_cols = F, gaps_col = 3)