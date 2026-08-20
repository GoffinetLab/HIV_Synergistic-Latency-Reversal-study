#### LOAD REQUIRED LIBRARIES ####------------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(knitr)
library(ggprism)
library(ggpubr)
library(rstatix)
library(gridExtra)
library(purrr)
library(magrittr)
library(patchwork)
library(tidyr)
library(readxl)

#### READ IN DATA ####-------------------------------------------------------------------------------------------------------

df.jleg_p24env = readxl::read_xlsx("JLEG.1 AZD+JQ1 synergy.xlsx")

#### PROCESS DATA FOR BLISS INDEPENDENCE TESTING ####------------------------------------------------------------------------

## perform calculation of bliss excess

for (param in c("p24_perc", "EnvGFP_perc", "surfaceEnv_perc")){
  
  df1 = df.jleg_p24env %>%
    dplyr::filter(Treatment %in% c("Mock", "AZD", "JQ1", "AZD+JQ1")) %>%
    dplyr::select(c(1:3, {{param}})) %>%
    rename(perc_positive={{param}}) %>%
    mutate(perc_positive = perc_positive / 100) %>%
    pivot_wider(names_from = Treatment, values_from = perc_positive) %>%
    mutate(across(AZD:`AZD+JQ1`, ~.x - Mock)) %>%
    dplyr::select(-Mock) %>%
    mutate(f_AJ_predicted = AZD + JQ1 - (AZD*JQ1))%>%
    mutate(bliss_excess = `AZD+JQ1` - f_AJ_predicted) %>%
    arrange(Clone) %>%
    mutate(param = {{param}})
  
  assign(paste0("df_bi.", param), value = df1)
  
  
}

df_bi.allparams = rbind(df_bi.p24_perc,
                        df_bi.EnvGFP_perc,
                        df_bi.surfaceEnv_perc)

#### PERFORM STATISTICS TEST FOR STAT SIGNIFICANCE OF BLISS EXCESS #### ---------------------------------

stat.test = df_bi.allparams %>%
  pivot_longer(cols = c("AZD+JQ1", "f_AJ_predicted"),
               names_to = "response_type",
               values_to = "response") %>%
  
  #ungroup()%>%
  group_by(Clone, param) %>%
  wilcox_test(response ~ response_type) %>%
  add_significance("p") %>%
  mutate(test = "wilcox")



stat.test_to_bind = stat.test %>%
  dplyr::select(Clone, param, p, p.signif)


df_bi.allparams.stats = plyr::join(df_bi.allparams, stat.test_to_bind)


#### PLOT RESULTS IN BAR GRAPH #### ----------------------------------------------------------------------

dat1 = df_bi.allparams.stats %>%
  mutate(param = factor(param, levels = c("p24_perc", "EnvGFP_perc", "surfaceEnv_perc")) )#%>%
#dplyr::filter(Clone =="JLEG2.1")

plot1 = ggplot(dat1, aes(param, bliss_excess))+
  stat_summary(fun = "mean", geom = "bar", colour="black", fill="#4657C7")+
  geom_point(aes(shape = as.factor(rep)), size=2,  position=position_jitter(height=0, width=0.1))+
  #scale_shape_manual(values=1:nlevels(as.factor(df$rep))) +
  scale_shape_manual(values=shape_set1) +
  stat_summary(fun.data = mean_sdl, geom = "errorbar", fun.args = list(mult = 1), width=0.15)+
  theme_prism()+
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust = 1))+
  scale_x_discrete(labels = c("p24CA", "Total EnvGFP", "Cell-surface EnvGFP"))+
  scale_y_continuous(limits = c(0, 0.8), expand = c(0,0),
                     guide = guide_prism_minor(),
                     breaks = scales::breaks_pretty(8))+
  geom_text(aes( y=0.65, label = p.signif, size=2))+
  theme(legend.position = "none") +
  facet_wrap(~Clone, labeller = labeller(Clone =
                                           c("JLEG2.1"="JLEG.1",
                                             "JLEG2.60"="JLEG.60",
                                             "JLEG2.74"="JLEG.74")))+
  theme(strip.text = element_text(size=14))+
  xlab("")+
  ylab("Bliss excess")
