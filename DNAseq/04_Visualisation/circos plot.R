# Genome version hg38, there are three circles in the picture. The outer circle represents the chromosome, the middle circle represents the positive strand gene, and the inner circle represents the negative strand gene:
# Use bar to display the relative position of the insertion point
# Use different colors of bars to represent gene biotypes
# Different directions represent strands (the upper circle of bars is the positive strand, and the lower circle of bars is the negative strand)
# The gene name hgnc symbol is marked with text

install.packages("circlize")
install.packages("ComplexHeatmap")
install.packages("openxlsx")
BiocManager::install("ComplexHeatmap")

library(tidyverse)
library(dplyr)
library(tidyr)
library(openxlsx)
library(circlize)
library(ComplexHeatmap)
library(grid)

# set digrid# set directory
setwd("~/R")

# read excel
data <- read.xlsx('input/01wgs_int.sites.xlsx', 1)%>%
    separate(chr.pos, c('chr', 'pos'), sep=':')%>% #sep chr and pos 
    mutate(chr = paste0('chr', chr))%>% #add chr
    mutate(pos = as.numeric(pos)) %>%
    mutate(hgnc_symbol = ifelse(is.na(hgnc_symbol), '', hgnc_symbol)) %>%
    dplyr::select(-clone)%>% #remove lib.no
    distinct

# set chromosome colours
colors <- c('#42150AFF', '#72190EFF', '#A5170EFF', '#DC050CFF', 
    '#E65518FF', '#E8601CFF', '#EE8026FF', '#F1932DFF', 
    '#F4A736FF', '#F6C141FF', '#F7CB45FF', '#F7F056FF', 
    '#CAE0ABFF', '#90C987FF', '#4EB265FF', '#7BAFDEFF', 
    '#5289C7FF', '#437DBFFF', '#1965B0FF', '#882E72FF', 
    '#994F88FF', '#AA6F9EFF', '#BA8DB4FF', '#CAACCBFF') %>%
    setNames(c(paste0('chr', 1:22), 'chrX', 'chrY'))

# set gene biotype colour
biotype_colors <- c('protein_coding' = 'black', 'lncRNA' = 'green', 'lncRNA+pseudogene' = 'orange')

# set intron/exon colour
intron <- c('Intron' = 'blue', 'Exon (UTR)' = '#F7F056FF', 'Exon (CDS)' = '#EE8026FF' )

# plot graph
#pdf(file = 'circos_10clones.pdf', width = 7, height = 7)

#jpeg(file = 'circos_10clones.jpg', width = 7, height = 7, units = "in", res = 300)

#svg(file = 'circos.svg', width = 8, height = 8, pointsize = 11)

#postscript("circos.eps", width = 8, height = 8, onefile = FALSE, paper = "special")

png("circos.png", width = 8, height = 8, units = "in", res = 300)


# set graphic parameters
par(mar = c(8, 4, 4, 2), xpd = NA)
circos.par(cell.padding = c(0, 0, 0, 0), start.degree = 90)

# set the graph into parts
# if only show chromosomes that are in the dataset  add this  index = unique(data$chr)
circos.initializeWithIdeogram(species = "hg38", 
    plotType = c("ideogram"))

# create plotting regions
# create inner chromosome track 
# circos.rect() give you dimensions of the plotting region and you can customize plot regions directly 
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    chr = CELL_META$sector.index
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    circos.rect(xlim[1], 0, xlim[2], 1, col = colors[chr], border = NA)
    circos.text(mean(xlim), mean(ylim), gsub('chr', '', chr), cex = 0.9, col = "white",
        facing = "downward", niceFacing = TRUE)
}, track.height = 0.15, bg.border = NA)

# add integration sites
# update plotting regions
# add tracks
# trial with lines
# 1st track - positive strand
data1 <- dplyr::filter(data, strand == '+')
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {}, track.height = 0.1, bg.border = NA)

circos.trackLines(
  data1$chr, data1$pos, rep(1, nrow(data1)),
  track.index = get.current.track.index(),
  col = "darkblue",
  lwd = 2, type = "h", baseline = "bottom"
)


# 2nd track - negative strand
data2 <- dplyr::filter(data, strand == '-')
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {}, track.height = 0.1, bg.border = NA)

circos.trackLines(
  data2$chr, data2$pos, rep(0, nrow(data2)),
  track.index = get.current.track.index(),
  col = "darkred",
  lwd = 2, type = "h", baseline = "top"
)

# show gene name

circos.genomicLabels(data[, c(3, 7, 8, 6)],  # Correct column selection
                     labels.column = 4,  # "hgnc_symbol" is now the 4th column in this subset
                     side = "inside",
                     facing = 'reverse.clockwise',
                     niceFacing = TRUE,
                     cex = 0.75, 
                     col = biotype_colors[data$gene_biotype])


# adding legend
lgd_biotype <- Legend(
  title = "Gene biotype",
  labels = c("Protein coding", "Long non-coding RNA", "Pseudogene / lncRNA"),
  legend_gp = gpar(fill = c("black", "green", "orange"))
)

lgd_strand <- Legend(
  title = "Strands",
  labels = c("Positive (Outer layer)", "Negative (Inner layer)"),
  type = "lines",
  legend_gp = gpar(col = c("darkred", "darkblue"), lwd = 2)
)

combined_lgd <- packLegend(lgd_biotype, lgd_strand, direction = "horizontal")

draw(
  combined_lgd,
  x = unit(0.5, "npc"),
  y = unit(0.1, "npc"),
  just = c("center", "bottom")
)

dev.off()
circos.clear()

