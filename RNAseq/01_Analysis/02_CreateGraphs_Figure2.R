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


## create helper functions

ggsave_png.svg = function(plot.obj, h=5, w =5, folder, name){
  
  ## save png
  ggsave(plot = plot.obj, filename = paste0(folder, "/", name, ".png"),
         height = h, width = w)
  ## save svg
  ggsave(plot = plot.obj, filename = paste0(folder, "/", name, ".svg"),
         height = h, width = w)
  
}


save_rds.xlsx = function(df, name, filepath){
  
  ## save as xlsx
  writexl::write_xlsx(df, path = paste0(filepath, "/", name, ".xlsx"))
  ## save as rds
  saveRDS(object = df, file = paste0(filepath, "/", name, ".rds"))
  
  
}


## set colours

cols_set = c("mock" = "#8F8F8F",
             "bryo" = "#E5625E",
             "jq1" = "#0C98E9",
             "bryo.jq1" = "#9B4CAC")

## load in custom list of gene panels

panels = readxl::read_excel("path/to/gene_panels.xlsx")


##### Figure 2 #################

##### Figure 2A PCA Plot ####

## read in dds
dds = readRDS(file="path/to/dds_clone_treat_group.rds")

## plot PCA
vst_dds <- vst(dds)

data <- plotPCA(vst_dds, intgroup = c( "group"), returnData=TRUE)
data$clone = gsub("(.*)_(.*)", replacement = "\\1", x = data$group)
data$treatment = gsub("(.*)_(.*)", replacement = "\\2", x = data$group)
data

percentVar <- round(100 * attr(data, "percentVar"))

pca <- ggplot(data, aes(PC1, PC2,  shape = clone, fill = treatment)) + 
  geom_point(size=3, stroke=1) +
  scale_fill_manual(values = c("mock"= "grey",
                               "bryo" = "#E5625E",
                               "jq1" = "#0C98E9",
                               "bryo.jq1" = "#9B4CAC"))+
  scale_shape_manual(values= c(24,21))+
  
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", face = "bold", size = 13.3))+
  theme(axis.text.y = element_text(colour = "black", face = "bold", size = 13.3))+
  theme(axis.title = element_text(colour = "black", face = "bold", size = 13.3))+
  geom_hline(linetype=2, linewidth = 0.75, yintercept = 0)+
  geom_vline(linetype=2, linewidth = 0.75, xintercept = 0)


#### Figure 2B ####

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

#### Create violin plots of up- and down-regulated genes

#### UNINFECTED 

dat2 = res_mega %>%
  dplyr::filter(padj < 0.05) %>%
  mutate(updown = ifelse(log2FoldChange <0, "down", "up")) %>%
  dplyr::filter(Clone == "Uninfect") %>%
  mutate(Treat = factor(Treat, levels =c("Bryo", "JQ1", "BryoJQ1")))



vln_un = ggplot(dat2, aes(fill = updown, y = log2FoldChange, x = Treat))+
  geom_violin(position=position_dodge(0), alpha=0.5, width=2) +
  geom_boxplot(aes(colour=updown), position=position_dodge(0), width=0.2, outlier.alpha = 0.5)+
  geom_boxplot(position=position_dodge(0), width=0.2, outlier.colour = "transparent")+
  scale_colour_manual(values= c(down = "blue", up="darkorange"))+
  scale_fill_manual(values= c(down = "blue", up="darkorange"))+
  theme_bw()+
  ylim(-10,10)+
  theme(axis.text.x = element_text(colour = "black", face = "bold", size = 13.3, angle=45, hjust=1, vjust = 1),
        axis.text.y = element_text(colour = "black", face = "bold", size = 13.3),
        axis.title.x = element_blank(),
        axis.title.y = element_text(colour = "black", face = "bold", size=13.3),
        #axis.title.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))+
  scale_x_discrete(labels = c("Bryo vs Mock", "JQ1 vs Mock", "Bryo+JQ1 vs Mock"))+
  geom_hline(yintercept = 0, linetype=2)+
  ggtitle("Uninfected")


##### JLEG2.1 

dat1 = res_mega %>%
  dplyr::filter(padj < 0.05) %>%
  mutate(updown = ifelse(log2FoldChange <0, "down", "up")) %>%
  dplyr::filter(Clone == "JLEG2.1") %>%
  mutate(Treat = factor(Treat, levels =c("Bryo", "JQ1", "BryoJQ1")))



vln_j2.1 = ggplot(dat1, aes(fill = updown, y = log2FoldChange, x = Treat))+
  geom_violin(position=position_dodge(0), alpha=0.5, width=2) +
  geom_boxplot(aes(colour=updown), position=position_dodge(0), width=0.2, outlier.alpha = 0.5)+
  geom_boxplot(position=position_dodge(0), width=0.2, outlier.colour = "transparent")+
  scale_colour_manual(values= c(down = "blue", up="darkorange"))+
  scale_fill_manual(values= c(down = "blue", up="darkorange"))+
  theme_bw()+
  ylim(-10,10)+
  theme(axis.text.x = element_text(colour = "black", face = "bold", size = 13.3, angle=45, hjust=1, vjust = 1),
        #axis.text.y = element_text(colour = "black", face = "bold", size = 10),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        #axis.title.y = element_text(colour = "black", face = "bold", size=12),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))+
  scale_x_discrete(labels = c("Bryo vs Mock", "JQ1 vs Mock", "Bryo+JQ1 vs Mock"))+
  geom_hline(yintercept = 0, linetype=2)+
  ggtitle("JLEG2.1")


#### Figure 2C ####

## create Venn diagrams for up and downregulated genes for JLEG2.1 clone

#### Upregulated genes

## JLEG2.1
bryo = deg_list$jleg2.1.bryo_v_mock %>%
  dplyr::filter(padj < 0.01 & log2FoldChange > 0) %>%
  pull(gene)
jq1 = deg_list$jleg2.1.jq1_v_mock %>%
  dplyr::filter(padj < 0.01 & log2FoldChange > 0) %>%
  pull(gene)
bryo.jq1 = deg_list$jleg2.1.bryo.jq1_v_mock %>%
  dplyr::filter(padj < 0.01 & log2FoldChange > 0) %>%
  pull(gene)

overlap_genes = list("Bryo" = bryo,
                     "JQ1" = jq1,
                     "Bryo+JQ1" = bryo.jq1)

## create venn diagram

plot.venn_j2.1 = ggvenn(stroke_alpha = 1,fill_alpha = 0.6,
                        overlap_genes, 
                        fill_color = c("#E5625E", "#0C98E9", "#9B4CAC"),
                        stroke_size = 0.5, set_name_size = 6,
                        text_size = 4.5,
                        show_percentage = F
)

#### Downregulated genes

## JLEG2.1
bryo = deg_list$jleg2.1.bryo_v_mock %>%
  dplyr::filter(padj < 0.01 & log2FoldChange < 0) %>%
  pull(gene)
jq1 = deg_list$jleg2.1.jq1_v_mock %>%
  dplyr::filter(padj < 0.01 & log2FoldChange < 0) %>%
  pull(gene)
bryo.jq1 = deg_list$jleg2.1.bryo.jq1_v_mock %>%
  dplyr::filter(padj < 0.01 & log2FoldChange < 0) %>%
  pull(gene)

overlap_genes = list("Bryo" = bryo,
                     "JQ1" = jq1,
                     "Bryo+JQ1" = bryo.jq1)

## create venn diagram

plot.venn_j2.1 = ggvenn(stroke_alpha = 1,fill_alpha = 0.6,
                        overlap_genes, 
                        fill_color = c("#E5625E", "#0C98E9", "#9B4CAC"),
                        stroke_size = 0.5, set_name_size = 6,
                        text_size = 4.5,
                        show_percentage = F
                        
)


#### Figure 2D ####

## create Venn diagrams for up and downregulated genes showing overlap for JLEG2.1 and uninfected clones for each trteatment individually

#### Upregulated genes

#BRYO
un.bryo_up = deg_list$jleg2.Un1.bryo_v_mock %>%
  dplyr::filter(padj < 0.01 & log2FoldChange > 0) %>%
  pull(gene)
j2.1_bryo_up = deg_list$jleg2.1.bryo_v_mock %>%
  dplyr::filter(padj < 0.01 & log2FoldChange > 0) %>%
  pull(gene)

overlap_genes = list("JLEG.1" = un.bryo_up,
                     "Uninfected" = j2.1_bryo_up)

plot.venn_clonalint_bryo_up = ggvenn(stroke_alpha = 0,fill_alpha = 0.6,
                                     overlap_genes, 
                                     fill_color = c("#E5625E", "white"),
                                     stroke_size = 0.5, set_name_size = 6, show_percentage=F, text_size=5)
plot.venn_clonalint_bryo_up

#JQ1
un.jq1_up = deg_list$jleg2.Un1.jq1_v_mock %>%
  dplyr::filter(padj < 0.01 & log2FoldChange > 0) %>%
  pull(gene)
j2.1_jq1_up = deg_list$jleg2.1.jq1_v_mock %>%
  dplyr::filter(padj < 0.01 & log2FoldChange > 0) %>%
  pull(gene)

overlap_genes = list("JLEG.1" = un.jq1_up,
                     "Uninfected" = j2.1_jq1_up)

plot.venn_clonalint_jq1_up = ggvenn(stroke_alpha = 0,fill_alpha = 0.6,
                                    overlap_genes, 
                                    fill_color = c("#0C98E9", "white"),
                                    stroke_size = 0.5, set_name_size = 6, show_percentage=F, text_size=5)
plot.venn_clonalint_jq1_up

#BRYO+JQ1
un.bryo.jq1_up = deg_list$jleg2.Un1.bryo.jq1_v_mock %>%
  dplyr::filter(padj < 0.01 & log2FoldChange > 0) %>%
  pull(gene)
j2.1_bryo.jq1_up = deg_list$jleg2.1.bryo.jq1_v_mock %>%
  dplyr::filter(padj < 0.01 & log2FoldChange > 0) %>%
  pull(gene)

overlap_genes = list("JLEG.1" = un.bryo.jq1_up,
                     "Uninfected" = j2.1_bryo.jq1_up)

plot.venn_clonalint_bryo.jq1_up = ggvenn(stroke_alpha = 0,fill_alpha = 0.6,
                                         overlap_genes, 
                                         fill_color = c("#9B4CAC", "white"),
                                         stroke_size = 0.5, set_name_size = 6, show_percentage=F, text_size=5)
plot.venn_clonalint_bryo.jq1_up


#### Downregulated genes

un.bryo_down = deg_list$jleg2.Un1.bryo_v_mock %>%
  dplyr::filter(padj < 0.01 & log2FoldChange < 0) %>%
  pull(gene)
j2.1_bryo_down = deg_list$jleg2.1.bryo_v_mock %>%
  dplyr::filter(padj < 0.01 & log2FoldChange < 0) %>%
  pull(gene)

overlap_genes = list("JLEG.1" = un.bryo_down,
                     "Uninfected" = j2.1_bryo_down)

plot.venn_clonalint_bryo_down = ggvenn(stroke_alpha = 0,fill_alpha = 0.6,
                                       overlap_genes, 
                                       fill_color = c("#E5625E", "white"),
                                       stroke_size = 0.5, set_name_size = 6, show_percentage=F, text_size=5)
plot.venn_clonalint_bryo_down

#JQ1
un.jq1_down = deg_list$jleg2.Un1.jq1_v_mock %>%
  dplyr::filter(padj < 0.01 & log2FoldChange < 0) %>%
  pull(gene)
j2.1_jq1_down = deg_list$jleg2.1.jq1_v_mock %>%
  dplyr::filter(padj < 0.01 & log2FoldChange < 0) %>%
  pull(gene)

overlap_genes = list("JLEG.1" = un.jq1_down,
                     "Uninfected" = j2.1_jq1_down)

plot.venn_clonalint_jq1_down = ggvenn(stroke_alpha = 0,fill_alpha = 0.6,
                                      overlap_genes, 
                                      fill_color = c("#0C98E9", "white"),
                                      stroke_size = 0.5, set_name_size = 6, show_percentage=F, text_size=5)
plot.venn_clonalint_jq1_down

#BRYO+JQ1
un.bryo.jq1_down = deg_list$jleg2.Un1.bryo.jq1_v_mock %>%
  dplyr::filter(padj < 0.01 & log2FoldChange < 0) %>%
  pull(gene)
j2.1_bryo.jq1_down = deg_list$jleg2.1.bryo.jq1_v_mock %>%
  dplyr::filter(padj < 0.01 & log2FoldChange < 0) %>%
  pull(gene)

overlap_genes = list("JLEG.1" = un.bryo.jq1_down,
                     "Uninfected" = j2.1_bryo.jq1_down)

plot.venn_clonalint_bryo.jq1_down = ggvenn(stroke_alpha = 0,fill_alpha = 0.6,
                                           overlap_genes, 
                                           fill_color = c("#9B4CAC", "white"),
                                           stroke_size = 0.5, set_name_size = 6, show_percentage=F, text_size=5)
plot.venn_clonalint_bryo.jq1_down


#### Figure 2E ####

## load previously created RDS object of gene set enrichment results for reactome database
gse_reactome.groupmodel = readRDS("path/to/gse_reactome.groupmodel.rds")

## convert results to single dataframe for graphing
pway.mega = gse_reactome.groupmodel %>%
  mutate(stat = ifelse(p.adjust < 0.05, "Significant", "NS"))
contrasts = unique(pway.mega$group)

contrasts_nointeract = contrasts[!grepl(pattern = "interact.bryo_jq1", contrasts)]
contrasts_nointeract

## define function to plot gene set enrichment results

## plot pathways on one dotplot
pway_comboDot <- function(contrasts_use = c(unique(pway.mega$group)),
                          num_topPways = 20,
                          minimum_NES = 1,
                          max_padj = 0.001,
                          name = "Pways_ComboDot_AllContrasts",
                          width=8,
                          height=10.5){
  
  ### PWAYS UP # select top n pathways by NES
  pways_up <- pway.mega %>%
    dplyr::filter(group %in% contrasts_use) %>%
    dplyr::filter(p.adjust < max_padj, 
                  NES > minimum_NES,
                  stat == "Significant") %>%
    group_by(group) %>%
    arrange(desc(NES)) %>%
    slice_head(n = num_topPways) %>%
    pull(var= Description)
  
  pways_up <- unique(pways_up) ## get names of pathways to plot
  
  
  
  ## cluster dataframe for better graphical representation
  mat <- pway.mega %>%
    dplyr::filter(group %in% contrasts_use) %>%
    subset(Description %in% pways_up) %>%
    dplyr::select(Description, NES, group) %>%
    pivot_wider(names_from = group, values_from = NES) %>%
    column_to_rownames(var = "Description")
  
  ## scale and cluster
  mat <- scale(mat)
  dist_mat <- dist(mat, method = 'euclidean')
  hclust_avg <- hclust(dist_mat, method = 'average')
  ## get plotting order
  label_order <- rownames(mat)[hclust_avg$order] 
  
  ## plot pways UP
  plot_up_df <- pway.mega %>%
    dplyr::filter(group %in% contrasts_use, stat == "Significant") %>%
    mutate(group = factor(group, levels = contrasts_use)) %>%
    subset(Description %in% pways_up) %>%
    mutate(Description = factor(Description, levels =label_order ))#%>%
  #arrange(match(rownames(.), label_order))# %>%
  
  plot_up <- ggplot(plot_up_df, aes(group, Description, color=factor(stat, levels = c("NS", "Significant")))) +
    geom_point(aes(size=GeneRatio, fill= NES, stroke=1), shape=21) +
    scale_fill_gradient(low='#0000FF',high='#FFFF00')+
    theme_bw() +
    # theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    scale_color_manual(values=c("black", "black"), name="Significance") +
    scale_x_discrete(drop=F)+
    #ggtitle("Pathways UP")+
    theme(axis.text.x = element_text(colour = "black", face = "bold", size = 13.3),
          axis.text.y = element_text(colour = "black", face = "bold", size = 10))
  ## save as PDF and svg
  ggsave(filename = paste0("path/to/dotplot_gse.reactome_", name, "_UP.svg"),
         plot = plot_up,
         width=width, height=height)
  
  
  ## PWAYS DOWN
  pways_down <- pway.mega %>%
    dplyr::filter(group %in% contrasts_use) %>%
    dplyr::filter(p.adjust < max_padj, 
                  NES < -minimum_NES,
                  stat == "Significant") %>%
    group_by(group) %>%
    arrange(NES) %>%
    slice_head(n = num_topPways) %>%
    pull(var= Description)
  
  pways_down <- unique(pways_down) ## get names of pathways to plot
  
  ## cluster dataframe for better graphical representation
  mat <- pway.mega %>%
    dplyr::filter(group %in% contrasts_use) %>%
    subset(Description %in% pways_down) %>%
    dplyr::select(Description, NES, group) %>%
    pivot_wider(names_from = group, values_from = NES) %>%
    column_to_rownames(var = "Description")
  
  ## scale and cluster
  mat <- scale(mat)
  dist_mat <- dist(mat, method = 'euclidean')
  dist_mat<- as.dist(dist_mat)
  ## some NAs created
  dist_mat[is.na(dist_mat)] <- 0
  dist_mat[is.nan(dist_mat)] <- 0
  sum(is.infinite(dist_mat))  # THIS SHOULD BE 0
  hclust_avg <- hclust(dist_mat, method = 'average')
  
  ## get plotting order
  label_order <- rownames(mat)[hclust_avg$order] 
  
  ## plot pways down
  plot_down_df <- pway.mega %>%
    dplyr::filter(group %in% contrasts_use, stat == "Significant") %>%
    mutate(group = factor(group, levels = contrasts_use)) %>%
    subset(Description %in% pways_down) %>%
    mutate(Description = factor(Description, levels =label_order ))#%>%
  #arrange(match(rownames(.), label_order))# %>%
  
  plot_down <- ggplot(plot_down_df, aes(group, Description, color=factor(stat, levels = c("NS", "Significant")))) +
    geom_point(aes(size=GeneRatio, fill= NES, stroke=1), shape=21) +
    scale_fill_gradient(low='#0000FF',high='#FFFF00')+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    scale_color_manual(values=c("black", "black"), name="Significance") +
    ggtitle("Pathways DOWN")+
    theme(axis.text.x = element_text(colour = "black", face = "bold", size = 13.3),
          axis.text.y = element_text(colour = "black", face = "bold", size = 10))
  
  ## save as PDF and svg
  ggsave(filename = paste0("path/to/dotplot_gse.reactome_", name, "_DOWN.svg"),
         plot = plot_down,
         width=width, height=height)
  
}

## creat dotplot with 10 highest scoring up and downregulated pathways

pway_comboDot(, name = "allcontrasts", num_topPways = 10, height =9, width =9)




