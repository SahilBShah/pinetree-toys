library(cowplot)
library(dplyr)
library(gridExtra)

setwd("/home/wilkelab/pinetree-toys/three_genes_evolution/")

original_data <- read.table("three_genes_test_file.tsv", header=TRUE)
original_data <- filter(original_data, species == "proteinX" | species == "proteinY" | species == "proteinZ")
original_line_plot <-  ggplot(original_data, aes(fill=species, color=species, x=time, y=transcript)) + geom_line(stat="identity") + ggtitle("Original Genome Function")

gen0_data <- read.table("gen_0_data.tsv", header=TRUE)
gen0_data <- filter(gen0_data, species == "proteinX" | species == "proteinY" | species == "proteinZ")
gen0_line_plot <-  ggplot(gen0_data, aes(fill=species, color=species, x=time, y=transcript)) + geom_line(stat="identity") + ggtitle("Genome Function at Generation 0")
  
gen500_data <- read.table("gen_500_data.tsv", header=TRUE)
gen500_data <- filter(gen500_data, species == "proteinX" | species == "proteinY" | species == "proteinZ")
gen500_line_plot <-  ggplot(gen500_data, aes(fill=species, color=species, x=time, y=transcript)) + geom_line(stat="identity") + ggtitle("Genome Function at Generation 500")

gen1000_data <- read.table("gen_1000_data.tsv", header=TRUE)
gen1000_data <- filter(gen500_data, species == "proteinX" | species == "proteinY" | species == "proteinZ")
gen1000_line_plot <-  ggplot(gen500_data, aes(fill=species, color=species, x=time, y=transcript)) + geom_line(stat="identity") + ggtitle("Genome Function at Generation 1000")

gen_final1_data <- read.table("best_replicated_gen_446.tsv", header=TRUE)
gen_final1_data <- filter(gen_final1_data, species == "proteinX" | species == "proteinY" | species == "proteinZ")
gen_final1_line_plot <-  ggplot(gen_final1_data, aes(fill=species, color=species, x=time, y=transcript)) + geom_line(stat="identity") + ggtitle("Genome Function of Best Genome Found")

gen_final2_data <- read.table("best_replicated_gen_2.tsv", header=TRUE)
gen_final2_data <- filter(gen_final2_data, species == "proteinX" | species == "proteinY" | species == "proteinZ")
gen_final2_line_plot <-  ggplot(gen_final2_data, aes(fill=species, color=species, x=time, y=transcript)) + geom_line(stat="identity") + ggtitle("Genome Function of Best Genome 2 Found")
  
grid.arrange(original_line_plot, gen0_line_plot, gen500_line_plot, gen1000_line_plot, gen_final1_line_plot, nrow=5)
grid.arrange(original_line_plot, gen_final1_line_plot, gen_final2_line_plot)

plot.new()
