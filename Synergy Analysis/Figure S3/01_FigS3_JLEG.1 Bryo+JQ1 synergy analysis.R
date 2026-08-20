#### LOAD REQUIRED LIBRARIES ####------------------------------------------------------------------------

library(ggplot2)
library(knitr)
library(ggprism)
library(ggpubr)
library(dplyr)
library(BIGL)

#### READ IN DATA ####------------------------------------------------------------------------------------

df.syn_all = readxl::read_xlsx("JLEG.1 Bryo+JQ1 synergy.xlsx")

#### SIMULATE NULL DISTRIBUTION

df.syn_all.mean = df.syn_all %>%
  dplyr::filter(Clone == "JLEG2.1" & Treatment == "Bryo_JQ1") %>%
  #mutate(dil_combo = paste0(Bryo_dil,"_", JQ1_dil)) %>%
  group_by(Clone, Treatment, Bryo_conc, JQ1_conc, Bryo_dil, JQ1_dil) %>%
  summarise(mean.p24_perc = mean(p24_perc),
            mean.EnvGFP_perc = mean(EnvGFP_perc),
            mean.surfaceEnv_perc = mean(surfaceEnv_perc)) 


## simulate null distribution by matrix addition

for (param in c("mean.p24_perc", "mean.EnvGFP_perc", "mean.surfaceEnv_perc")){
  
  print(param)
  
  
  ## create empty list to write to 
  append_nullres <- c()
  temp <- df.syn_all.mean %>%
    filter(Treatment == "Bryo_JQ1")
  ## pull Bryo results for JQ1=0
  bryoNull <- temp %>%
    filter(JQ1_dil == 0) %>%
    pull(param)
  ## pull JQ1 results for Bryo=0
  jq1Null <- temp %>%
    filter(Bryo_dil == 0) %>%
    pull(param)
  
  # start counter at 0
  j <- 0
  addNull_res <- list()
  # matrix addition = add list of JQ1 results to result at each Bryo conc
  for (i in bryoNull){
    j <- j +1
    addNull <- i + jq1Null
    addNull_res[[j]] <- addNull
    
  }
  ## unlist 
  append_nullres <- append(x= append_nullres, unlist(addNull_res))
  ## add to df
  ## name for new column
  null_name = paste0(param, ".null")
  df.syn_all.mean[null_name] = append_nullres
  
  
}


#### PLOT SYNERGY LANDSCAPES ####------------------------------------------------------------------------------------

#### p24 #######

p24.1 = df.syn_all.mean %>%
  
  dplyr::select(1:6, mean.p24_perc) %>%
  mutate(distrib = "observed") %>%
  rename(reading=mean.p24_perc) %>%
  mutate(param = "p24_perc")

p24.2 = df.syn_all.mean %>%
  
  dplyr::select(1:6, mean.p24_perc.null) %>%
  rename(reading=mean.p24_perc.null) %>%
  mutate(distrib = "null") %>%
  mutate(param = "p24_perc") %>%
  rbind(., p24.1)

#### EnvGFP #######

EnvGFP.1 = df.syn_all.mean %>%
  
  dplyr::select(1:6, mean.EnvGFP_perc) %>%
  mutate(distrib = "observed") %>%
  rename(reading=mean.EnvGFP_perc) %>%
  mutate(param = "EnvGFP_perc")

EnvGFP.2 = df.syn_all.mean %>%
  
  dplyr::select(1:6, mean.EnvGFP_perc.null) %>%
  rename(reading=mean.EnvGFP_perc.null) %>%
  mutate(distrib = "null") %>%
  mutate(param = "EnvGFP_perc") %>%
  rbind(., EnvGFP.1)

#### surfaceEnv #######

surfaceEnv.1 = df.syn_all.mean %>%
  
  dplyr::select(1:6, mean.surfaceEnv_perc) %>%
  mutate(distrib = "observed") %>%
  rename(reading=mean.surfaceEnv_perc) %>%
  mutate(param = "surfaceEnv_perc")

surfaceEnv.2 = df.syn_all.mean %>%
  
  dplyr::select(1:6, mean.surfaceEnv_perc.null) %>%
  rename(reading=mean.surfaceEnv_perc.null) %>%
  mutate(distrib = "null") %>%
  mutate(param = "surfaceEnv_perc") %>%
  rbind(., surfaceEnv.1)

### bind all ####

df.meanparam = rbind(p24.2,
                     EnvGFP.2,
                     surfaceEnv.2) %>%
  mutate(param = factor(param, levels = c("p24_perc", "EnvGFP_perc", "surfaceEnv_perc")))


## save as raster


p2 <- ggplot(df.meanparam, aes(Bryo_dil, JQ1_dil,  z=reading)) +
  geom_contour_filled(bins=10, alpha=0.7)+
  geom_contour(colour = "black", bins=10)+
  scale_fill_brewer(type = 'div', palette = 'Spectral', direction = -1)+
  facet_wrap(~distrib+ param)+
  theme_bw()+
  scale_x_continuous(breaks=c(seq(0, 7, by=1)))+
  scale_y_continuous(breaks=c(seq(0, 7, by=1)))

p2

#### PERFORM SYNERGY ANALYSIS WITH BIGL PACKAGE #### ------------------------------------------------------

# p24 reading synergy____________________________________________________________________________

## marginal fit
df.syntest_p24 <- df.syn_all %>%
  mutate(Bryo_conc = Bryo_conc * 1000, JQ1_conc = JQ1_conc * 1000) %>%
  filter(Clone == "JLEG2.1", Treatment == "Bryo_JQ1") %>%
  ungroup() %>%
  select(Bryo_conc, JQ1_conc, p24_perc) %>%
  rename("d1" = "Bryo_conc") %>%
  rename("d2" = "JQ1_conc") %>%
  rename("effect" = "p24_perc")

margFit_p24 <- fitMarginals(df.syntest_p24, method = "optim", names = c("Bryo", "JQ1"))

## Bliss Independence test
rsb_p24 <- fitSurface(df.syntest_p24, margFit_p24, 
                      null_model = "bliss",
                      B.CP = 50, statistic = "maxR", parallel = FALSE)

## synergy contour plot
p24.syn= summary(rsb_p24)$maxR[[2]] %>%
  rename("Bryo_conc"="d1", "JQ1_conc"="d2") %>%
  # mutate(treat_combo = paste0(Bryo_conc, "_", JQ1_conc))
  mutate(JQ1_dil = ifelse(JQ1_conc == 80, 3,
                          ifelse(JQ1_conc == 400, 4,
                                 ifelse(JQ1_conc == 2000, 5,
                                        ifelse(JQ1_conc == 10000, 6,
                                               ifelse(JQ1_conc == 50000, 7,
                                                      ifelse(JQ1_conc == 16, 2,
                                                             ifelse(JQ1_conc == 3.2, 1,NA)))))))) %>%
  mutate(Bryo_dil = ifelse(Bryo_conc == 0.8, 3,
                           ifelse(Bryo_conc == 4, 4,
                                  ifelse(Bryo_conc == 20, 5,
                                         ifelse(Bryo_conc == 100, 6,
                                                ifelse(Bryo_conc == 500, 7, 
                                                       ifelse(Bryo_conc == 0.16, 2,
                                                              ifelse(Bryo_conc == 0.032, 1,1)))))))) %>%
  rename("pvalue"="p-value")

## draw 


## synergy
syn.hull_p001 = p24.syn %>%
  filter(call=="Syn" & pvalue < 0.001) %>%
  slice(chull(Bryo_dil, JQ1_dil)) %>%
  mutate(sig = "syn_001")
syn.hull_p005 = p24.syn %>%
  filter(call=="Syn" &  pvalue <0.05) %>%
  slice(chull(Bryo_dil, JQ1_dil))%>%
  mutate(sig = "syn_005")
syn.hull_p01 = p24.syn %>%
  filter(call=="Syn" &  pvalue <0.01) %>%
  slice(chull(Bryo_dil, JQ1_dil))%>%
  mutate(sig = "syn_01")

## antagonism
ant.hull_p001 = p24.syn %>%
  filter(call=="Ant" & pvalue < 0.001) %>%
  slice(chull(Bryo_dil, JQ1_dil)) %>%
  mutate(sig = "ant_001")
ant.hull_p005 = p24.syn %>%
  filter(call=="Ant" &  pvalue <0.05) %>%
  slice(chull(Bryo_dil, JQ1_dil))%>%
  mutate(sig = "ant_005")
ant.hull_p01 = p24.syn %>%
  filter(call=="Ant" &  pvalue <0.01) %>%
  slice(chull(Bryo_dil, JQ1_dil))%>%
  mutate(sig = "ant_01")


syn_ant = rbind(syn.hull_p001, syn.hull_p005, syn.hull_p01,
                ant.hull_p001, ant.hull_p005, ant.hull_p01) %>%
  mutate(sig = factor(sig, levels = c("ant_01", "ant_005", "ant_001",
                                      "syn_01", "syn_005", "syn_001")))


syn_plot.p24 =ggplot(p24.syn, aes(Bryo_dil, JQ1_dil)) +
  
  geom_polygon(data=syn_ant, alpha=0.6,
               aes(fill=sig))+
  scale_fill_manual(values = c("#F88379", "#DE3163","#8B0000",
                               "#ADD8E6", "#6495ED", "#0047AB"))+
  
  
  geom_point(aes(size=absR), alpha=0.8)+
  scale_size_continuous(limits=c(1,25))+
  theme_bw()+
  expand_limits(x=0, y=0)+
  scale_x_continuous( breaks=c(seq(0, 7, by=1)))+
  scale_y_continuous( breaks=c(seq(0, 7, by=1)))#+

# EnvGFP reading synergy ____________________________________________________________________________

df.syntest_EnvGFP <- df.syn_all %>%
  mutate(Bryo_conc = Bryo_conc * 1000, JQ1_conc = JQ1_conc * 1000) %>%
  filter(Clone == "JLEG2.1", Treatment == "Bryo_JQ1") %>%
  ungroup() %>%
  select(Bryo_conc, JQ1_conc, EnvGFP_perc) %>%
  rename("d1" = "Bryo_conc") %>%
  rename("d2" = "JQ1_conc") %>%
  rename("effect" = "EnvGFP_perc")

## marginal fit
margFit_EnvGFP <- fitMarginals(df.syntest_EnvGFP, method = "optim", names = c("Bryo", "JQ1"))

## Bliss Independence test
rsb_EnvGFP <- fitSurface(df.syntest_EnvGFP, margFit_EnvGFP, 
                         null_model = "bliss",
                         B.CP = 50, statistic = "maxR", parallel = FALSE)


## synergy contour
EnvGFP.syn= summary(rsb_EnvGFP)$maxR[[2]] %>%
  rename("Bryo_conc"="d1", "JQ1_conc"="d2") %>%
  # mutate(treat_combo = paste0(Bryo_conc, "_", JQ1_conc))
  mutate(JQ1_dil = ifelse(JQ1_conc == 80, 3,
                          ifelse(JQ1_conc == 400, 4,
                                 ifelse(JQ1_conc == 2000, 5,
                                        ifelse(JQ1_conc == 10000, 6,
                                               ifelse(JQ1_conc == 50000, 7,
                                                      ifelse(JQ1_conc == 16, 2,
                                                             ifelse(JQ1_conc == 3.2, 1,NA)))))))) %>%
  mutate(Bryo_dil = ifelse(Bryo_conc == 0.8, 3,
                           ifelse(Bryo_conc == 4, 4,
                                  ifelse(Bryo_conc == 20, 5,
                                         ifelse(Bryo_conc == 100, 6,
                                                ifelse(Bryo_conc == 500, 7, 
                                                       ifelse(Bryo_conc == 0.16, 2,
                                                              ifelse(Bryo_conc == 0.032, 1,1)))))))) %>%
  rename("pvalue"="p-value")

## draw convex hulls for synergy

## synergy
syn.hull_p001 = EnvGFP.syn %>%
  filter(call=="Syn" & pvalue < 0.001) %>%
  slice(chull(Bryo_dil, JQ1_dil)) %>%
  mutate(sig = "syn_001")
syn.hull_p005 = EnvGFP.syn %>%
  filter(call=="Syn" &  pvalue <0.05) %>%
  slice(chull(Bryo_dil, JQ1_dil))%>%
  mutate(sig = "syn_005")
syn.hull_p01 = EnvGFP.syn %>%
  filter(call=="Syn" &  pvalue <0.01) %>%
  slice(chull(Bryo_dil, JQ1_dil))%>%
  mutate(sig = "syn_01")

## antagonism
ant.hull_p001 = EnvGFP.syn %>%
  filter(call=="Ant" & pvalue < 0.001) %>%
  slice(chull(Bryo_dil, JQ1_dil)) %>%
  mutate(sig = "ant_001")
ant.hull_p005 = EnvGFP.syn %>%
  filter(call=="Ant" &  pvalue <0.05) %>%
  slice(chull(Bryo_dil, JQ1_dil))%>%
  mutate(sig = "ant_005")
ant.hull_p01 = EnvGFP.syn %>%
  filter(call=="Ant" &  pvalue <0.01) %>%
  slice(chull(Bryo_dil, JQ1_dil))%>%
  mutate(sig = "ant_01")


syn_ant = rbind(syn.hull_p001, syn.hull_p005, syn.hull_p01,
                ant.hull_p001, ant.hull_p005, ant.hull_p01) %>%
  mutate(sig = factor(sig, levels = c("ant_01", "ant_005", "ant_001",
                                      "syn_01", "syn_005", "syn_001")))


syn_plot.EnvGFP =ggplot(EnvGFP.syn, aes(Bryo_dil, JQ1_dil)) +
  
  geom_polygon(data=syn_ant, alpha=0.6,
               aes(fill=sig))+
  scale_fill_manual(values = c("#F88379", "#DE3163","#8B0000",
                               "#ADD8E6", "#6495ED", "#0047AB"))+
  
  
  geom_point(aes(size=absR), alpha=0.8)+
  scale_size_continuous(limits=c(1,25))+
  theme_bw()+
  expand_limits(x=0, y=0)+
  scale_x_continuous( breaks=c(seq(0, 7, by=1)))+
  scale_y_continuous( breaks=c(seq(0, 7, by=1)))

# Cell-surface EnvGFP reading synergy____________________________________________________________________________

df.syntest_surfaceEnv <- df.syn_all %>%
  mutate(Bryo_conc = Bryo_conc * 1000, JQ1_conc = JQ1_conc * 1000) %>%
  filter(Clone == "JLEG2.1", Treatment == "Bryo_JQ1") %>%
  ungroup() %>%
  select(Bryo_conc, JQ1_conc, surfaceEnv_perc) %>%
  rename("d1" = "Bryo_conc") %>%
  rename("d2" = "JQ1_conc") %>%
  rename("effect" = "surfaceEnv_perc")

## marginal fit
margFit_surfaceEnv <- fitMarginals(df.syntest_surfaceEnv, method = "optim", names = c("Bryo", "JQ1"))

## Bliss Independence test
rsb_surfaceEnv <- fitSurface(df.syntest_surfaceEnv, margFit_surfaceEnv, 
                             null_model = "bliss",
                             B.CP = 50, statistic = "maxR", parallel = FALSE)

## synergy contour
surfaceEnv.syn= summary(rsb_surfaceEnv)$maxR[[2]] %>%
  rename("Bryo_conc"="d1", "JQ1_conc"="d2") %>%
  # mutate(treat_combo = paste0(Bryo_conc, "_", JQ1_conc))
  mutate(JQ1_dil = ifelse(JQ1_conc == 80, 3,
                          ifelse(JQ1_conc == 400, 4,
                                 ifelse(JQ1_conc == 2000, 5,
                                        ifelse(JQ1_conc == 10000, 6,
                                               ifelse(JQ1_conc == 50000, 7,
                                                      ifelse(JQ1_conc == 16, 2,
                                                             ifelse(JQ1_conc == 3.2, 1,NA)))))))) %>%
  mutate(Bryo_dil = ifelse(Bryo_conc == 0.8, 3,
                           ifelse(Bryo_conc == 4, 4,
                                  ifelse(Bryo_conc == 20, 5,
                                         ifelse(Bryo_conc == 100, 6,
                                                ifelse(Bryo_conc == 500, 7, 
                                                       ifelse(Bryo_conc == 0.16, 2,
                                                              ifelse(Bryo_conc == 0.032, 1,1)))))))) %>%
  rename("pvalue"="p-value")

## draw convex hulls for synergy

## synergy
syn.hull_p001 = surfaceEnv.syn %>%
  filter(call=="Syn" & pvalue < 0.001) %>%
  slice(chull(Bryo_dil, JQ1_dil)) %>%
  mutate(sig = "syn_001")
syn.hull_p005 = surfaceEnv.syn %>%
  filter(call=="Syn" &  pvalue <0.05) %>%
  slice(chull(Bryo_dil, JQ1_dil))%>%
  mutate(sig = "syn_005")
syn.hull_p01 = surfaceEnv.syn %>%
  filter(call=="Syn" &  pvalue <0.01) %>%
  slice(chull(Bryo_dil, JQ1_dil))%>%
  mutate(sig = "syn_01")

## antagonism
ant.hull_p001 = surfaceEnv.syn %>%
  filter(call=="Ant" & pvalue < 0.001) %>%
  slice(chull(Bryo_dil, JQ1_dil)) %>%
  mutate(sig = "ant_001")
ant.hull_p005 = surfaceEnv.syn %>%
  filter(call=="Ant" &  pvalue <0.05) %>%
  slice(chull(Bryo_dil, JQ1_dil))%>%
  mutate(sig = "ant_005")
ant.hull_p01 = surfaceEnv.syn %>%
  filter(call=="Ant" &  pvalue <0.01) %>%
  slice(chull(Bryo_dil, JQ1_dil))%>%
  mutate(sig = "ant_01")



syn_ant = rbind(syn.hull_p001, syn.hull_p005, syn.hull_p01,
                ant.hull_p001, ant.hull_p005, ant.hull_p01) %>%
  mutate(sig = factor(sig, levels = c("ant_01", "ant_005", "ant_001",
                                      "syn_01", "syn_005", "syn_001")))


syn_plot.surfaceEnv =ggplot(surfaceEnv.syn, aes(Bryo_dil, JQ1_dil)) +
  
  geom_polygon(data=syn_ant, alpha=0.6,
               aes(fill=sig))+
  scale_fill_manual(values = c("#F88379", "#DE3163","#8B0000",
                               "#ADD8E6", "#6495ED", "#0047AB"))+
  
  
  geom_point(aes(size=absR), alpha=0.8)+
  scale_size_continuous(limits=c(1,25))+
  theme_bw()+
  scale_x_continuous( breaks=c(seq(0, 7, by=1)), limits = c(0,7)) +
  scale_y_continuous( breaks=c(seq(0, 7, by=1)), limits = c(0,7))

