library(cowplot)
library(dplyr)
library(gridExtra)

install.packages("gridExtra")

setwd("/home/wilkelab/pinetree-toys/three_genes_evolution/")

original_data <- read.table("three_genes_test_file.tsv", header=TRUE)
original_data <- filter(original_data, species == "proteinX" | species == "proteinY" | species == "proteinZ")
original_line_plot <-  ggplot(original_data, aes(fill=species, color=species, x=time, y=transcript)) + geom_point(stat="identity") + geom_line(stat="identity") + ggtitle("Original Genome Function")

gen0_data <- read.table("gen_0_data.tsv", header=TRUE)
gen0_data <- filter(gen0_data, species == "proteinX" | species == "proteinY" | species == "proteinZ")
gen0_line_plot <-  ggplot(gen0_data, aes(fill=species, color=species, x=time, y=transcript)) + geom_point(stat="identity") + geom_line(stat="identity") + ggtitle("Genome Function at Generation 0")

gen500_data <- read.table("gen_500_data.tsv", header=TRUE)
gen500_data <- filter(gen500_data, species == "proteinX" | species == "proteinY" | species == "proteinZ")
gen500_line_plot <-  ggplot(gen500_data, aes(fill=species, color=species, x=time, y=transcript)) + geom_point(stat="identity") + geom_line(stat="identity") + ggtitle("Genome Function at Generation 500")

gen_final_data <- read.table("best_three_genes_replicated.tsv", header=TRUE)
gen_final_data <- filter(gen_final_data, species == "proteinX" | species == "proteinY" | species == "proteinZ")
gen_final_line_plot <-  ggplot(gen_final_data, aes(fill=species, color=species, x=time, y=transcript)) + geom_point(stat="identity") + geom_line(stat="identity") + ggtitle("Genome Function of Best Genome Found")

grid.arrange(original_line_plot, gen0_line_plot, gen500_line_plot, gen_final_line_plot, nrow=4)
