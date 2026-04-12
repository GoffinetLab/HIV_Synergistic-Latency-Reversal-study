suppressMessages(library(DESeq2))
library(ggplot2)
library(tidyverse)
library(dplyr)
library(AnnotationDbi)
library(EnsDb.Hsapiens.v86)
library(biomaRt)
library(apeglm)
library(msigdbr)
library(GSVA)
library(clusterProfiler)
#library(clusterProfiler, lib.loc = "C:/Program Files/R/R-4.1.0/library")
library(org.Hs.eg.db)
library(ReactomePA)
library(freerange)
library(reactome.db)

#---------------------------------------------------------------------

## Read in DESeq2 dds result

dds = readRDS("path/to/dds_clone_treat_group.rds")
resultsNames(dds)

###### DIFFERENTIAL EXPRESSION #######

## define function to do all DEGs (known to man)

## setup function to process DESEQ2 contrast results in consistent way
## convert to list of dataframes, referred to as res_list
## add col for ensID 

res_list = list()

process_res = function(res_df, name){
  res1 = lfcShrink(dds, res = res_df, type = "ashr")
  res2 = as.data.frame(res1) 
  res2$ensID = rownames(res2)
  res2$gene <- mapIds(EnsDb.Hsapiens.v86, keys= res2$ensID, column="SYMBOL", keytype="GENEID", multiVals="first")
  res2$gene = ifelse(is.na(res2$gene), res2$ensID, res2$gene)
  res_list[[name]] <- res2
  assign("res_list", res_list, envir = .GlobalEnv)
  
}

do_deg_all = function(dds_obj){
  
  res_list = list()
  
  ###### JLEG2.1 #######
  
  ## 1 JLEG2.1 Bryo v Mock ####
  res = results(dds, contrast=c(1,0,0,-1,0,0,0,0,0,0))
  process_res(res_df=res, name="jleg2.1.bryo_v_mock")
  
  ## 2 JLEG2.1 JQ1 v Mock ####
  res = results(dds, contrast=c(0,0,1,-1,0,0,0,0,0,0))
  process_res(res_df=res, name="jleg2.1.jq1_v_mock")
  
  ## 3 JLEG2.1 Bryo+JQ1 v Mock ####
  res = results(dds, contrast=c(0,1,0,-1,0,0,0,0,0,0))
  process_res(res_df=res, name="jleg2.1.bryo.jq1_v_mock")
  
  ## 4 JLEG2.1 Interact Bryo+JQ1 ####
  res = results(dds, contrast=c(-1,1,-1,1,0,0,0,0,0,0))
  process_res(res_df=res, name="jleg2.1.interact.bryo_jq1")
  
  ###### JLEG2.Un1 #######
  
  ## 1 JLEG2.Un1 Bryo v Mock ####
  res = results(dds, contrast=c(0,0,0,0,1,0,0,-1,0,0))
  process_res(res_df=res, name="jleg2.Un1.bryo_v_mock")
  
  ## 2 JLEG2.Un1 JQ1 v Mock ####
  res = results(dds, contrast=c(0,0,0,0,0,0,1,-1,0,0))
  
  process_res(res_df=res, name="jleg2.Un1.jq1_v_mock")
  
  ## 3 JLEG2.Un1 Bryo+JQ1 v Mock ####
  res = results(dds, contrast=c(0,0,0,0,0,1,0,-1,0,0))
  process_res(res_df=res, name="jleg2.Un1.bryo.jq1_v_mock")
  
  ## 4 JLEG2.Un1 Interact Bryo+JQ1 ####
  res = results(dds, contrast=c(0,0,0,0,-1,1,-1,1,0,0))
  process_res(res_df=res, name="jleg2.Un1.interact.bryo_jq1")
  
  ###### CLONAL INTERACTION #######
  
  ## 1 CLONAL.INTERACT Bryo v Mock ####
  res = results(dds, contrast=c(1,0,0,-1,-1,0,0,1,0,0))
  process_res(res_df=res, name="clonal.interact.bryo_v_mock")
  
  ## 2 CLONAL.INTERACT JQ1 v Mock ####
  res = results(dds, contrast=c(0,0,1,-1,0,0,-1,1,0,0))
  process_res(res_df=res, name="clonal.interact.jq1_v_mock")
  
  ## 3 CLONAL.INTERACT JQ1 v Mock ####
  res = results(dds, contrast=c(0,1,0,-1,0,-1,0,1,0,0))
  process_res(res_df=res, name="clonal.interact.bryo.jq1_v_mock")
}

do_deg_all(dds)

##### save res_list #####


saveRDS(res_list, file = "path/to/res_list.rds")




#### GENESET ENRICHMENT ######



# Gene Ontology:Biological Process #######################################################

## define function to do GOBP enrichment

res_list$jleg2.1.bryo_v_mock

input = res_list2
contr = "jleg2.1.bryo_v_mock"


res_list2 = res_list[1:2]


run_go = function(input, name){
  
  go_res_df = data.frame()
  
  for (contr in names(res_list)){
    
    go_input <- input[[contr]] %>%
      #assign ENTREZID ID from symbol
      mutate(entrezID = mapIds(org.Hs.eg.db, keys=.$ensID, column="ENTREZID", keytype="ENSEMBL", multiVals="first")) %>%
      group_by(gene)%>%
      arrange(desc(log2FoldChange)) %>%
      dplyr::filter(row_number() == 1) %>% ## multimappers, keep entry w/ highest FC
      ungroup() %>%
      dplyr::filter(padj < 0.05) %>%
      dplyr::select(c(entrezID, log2FoldChange)) %>%
      na.omit () %>% # remove unmapped
      arrange(desc(log2FoldChange)) %>%
      deframe()
    
    
    ## run GO enrichment
    goenr <- gseGO(geneList     = go_input,
                   OrgDb        = org.Hs.eg.db,
                   ont          = "BP",
                   minGSSize    = 20,
                   maxGSSize    = 500,
                   pvalueCutoff = 1,
                   verbose      = T)
    
    ## get result
    goenr_tidy <-goenr@result 
    
    if (nrow(goenr_tidy) != 0){
      ## clean and determine GeneRatio
      gene_count <- goenr_tidy %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)
      merged_res <- left_join(goenr_tidy, gene_count, by = "ID") %>% 
        mutate(GeneRatio = count/setSize) %>%
        mutate(zscore = scale(NES)) %>%
        mutate(group = contr)
      
      ## add goRes to dataframe
      go_res_df = rbind(go_res_df, merged_res)
      
    }
    
    
  }
  
  ## export to global environment
  assign(x = paste0("gse_GOBP.", name), value = go_res_df, envir = .GlobalEnv)
  saveRDS(object = go_res_df, file = paste0("path/to/","gse_GOBP.",
                                            name, ".rds"))
}


run_go(input = res_list, name = "groupmodel")

run_go(input = res_list2, name = "test")


#### REACTOME #####################################################################

run_reactome = function(input, name){
  
  go_res_df = data.frame()
  
  for (contr in names(res_list)){
    
    go_input <- input[[contr]] %>%
      #assign ENTREZID ID from symbol
      mutate(entrezID = mapIds(org.Hs.eg.db, keys=.$ensID, column="ENTREZID", keytype="ENSEMBL", multiVals="first")) %>%
      group_by(gene)%>%
      arrange(desc(log2FoldChange)) %>%
      dplyr::filter(row_number() == 1) %>% ## multimappers, keep entry w/ highest FC
      ungroup() %>%
      dplyr::filter(padj < 0.05) %>%
      dplyr::select(c(entrezID, log2FoldChange)) %>%
      na.omit () %>% # remove unmapped
      arrange(desc(log2FoldChange)) %>%
      deframe()
    ## run GO enrichment
    reactenr <- gsePathway(geneList= go_input,
                           maxGSSize    = 500,
                           pvalueCutoff = 1,
                           pAdjustMethod = "BH",
                           verbose      = T)
    
    ## get result
    goenr_tidy <-reactenr@result 
    
    if (nrow(goenr_tidy) != 0){
      ## clean and determine GeneRatio
      gene_count <- goenr_tidy %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)
      merged_res <- left_join(goenr_tidy, gene_count, by = "ID") %>% 
        mutate(GeneRatio = count/setSize) %>%
        mutate(zscore = scale(NES)) %>%
        mutate(group = contr)
      
      ## add goRes to dataframe
      go_res_df = rbind(go_res_df, merged_res)
      
    }
    
    
  }
  
  ## export to global environment
  assign(x = paste0("gse_reactome.", name), value = go_res_df, envir = .GlobalEnv)
  saveRDS(object = go_res_df, file = paste0("path/to/","gse_reactome.",
                                            name, ".rds"))
}


run_reactome(input = res_list, name = "groupmodel")





keytypes()
## function to run for single DEG res
run_go_singledf = function(input, name, contrast){
  
  go_res_df = data.frame()
  
  go_input <- input %>%
    #assign ENTREZID ID from symbol
    mutate(ensembl_id = mapIds(org.Hs.eg.db, keys=.$gene, column="ENTREZID", keytype="SYMBOL", multiVals="first")) %>%
    select(c(ensembl_id, avg_log2FC)) %>%
    na.omit () %>% # remove unmapped
    arrange(desc(avg_log2FC)) %>%
    deframe()
  ## run GO enrichment
  goenr <- gseGO(geneList     = go_input,
                 OrgDb        = org.Hs.eg.db,
                 ont          = "BP",
                 minGSSize    = 20,
                 maxGSSize    = 500,
                 pvalueCutoff = 1,
                 verbose      = T)
  
  ## get result
  goenr_tidy <-goenr@result 
  
  gene_count <- goenr_tidy %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)
  merged_res <- left_join(goenr_tidy, gene_count, by = "ID") %>% 
    mutate(GeneRatio = count/setSize) %>%
    mutate(zscore = scale(NES)) %>%
    mutate(group = contrast)
  
  
  ## add goRes to dataframe
  go_res_df = rbind(go_res_df, merged_res)
  
  
  ## export to global environment
  assign(x = paste0("gse_GOBP.", name), value = go_res_df, envir = .GlobalEnv)
  saveRDS(object = go_res_df, file = paste0("path/to/","gse_reactome.",
                                            name, ".rds"))
}

run_reactome(input = res_list, name = "groupmodel")

