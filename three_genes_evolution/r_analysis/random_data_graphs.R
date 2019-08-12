library(cowplot)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(gtools)

setwd("/home/sahil/pinetree-toys/three_genes_evolution/fixed_memory_leak/")





original_data <- read.table("random_data.tsv", header=TRUE)
original_data <- filter(original_data, species == "proteinX" | species == "proteinY" | species == "proteinZ")
original_line_plot <-  ggplot(original_data, aes(fill=species, color=species, x=time, y=transcript)) + geom_line(stat="identity") + ggtitle("Original Genome Function") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

gen_plots <- grid.arrange(original_line_plot, nrow=1, ncol=1)
ggsave("/home/sahil/pinetree-toys/three_genes_evolution/fixed_memory_leak/random_data_graph.jpg", gen_plots, width=20, height=10)
