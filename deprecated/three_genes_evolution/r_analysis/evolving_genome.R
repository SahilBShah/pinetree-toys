library(cowplot)
library(dplyr)
library(gridExtra)
library(ggplot2)

setwd("/home/wilkelab/pinetree-toys/three_genes_evolution/")

best1 <- read.table("best_replicated_gen_15.tsv", header=TRUE)
best1 <- filter(best1, species == "proteinX" | species == "proteinY" | species == "proteinZ")
best1_plot <-  ggplot(best1, aes(fill=species, color=species, x=time, y=transcript)) + geom_line(stat="identity") + ggtitle("Genome Function at First Evolution")

best2 <- read.table("best_replicated_gen_185.tsv", header=TRUE)
best2 <- filter(best2, species == "proteinX" | species == "proteinY" | species == "proteinZ")
best2_plot <-  ggplot(best2, aes(fill=species, color=species, x=time, y=transcript)) + geom_line(stat="identity") + ggtitle("Genome Function at Second Evolution")

best3 <- read.table("best_replicated_gen_186.tsv", header=TRUE)
best3 <- filter(best3, species == "proteinX" | species == "proteinY" | species == "proteinZ")
best3_plot <-  ggplot(best3, aes(fill=species, color=species, x=time, y=transcript)) + geom_line(stat="identity") + ggtitle("Genome Function at Third Evolution")

best4 <- read.table("best_replicated_gen_187.tsv", header=TRUE)
best4 <- filter(best4, species == "proteinX" | species == "proteinY" | species == "proteinZ")
best4_plot <-  ggplot(best4, aes(fill=species, color=species, x=time, y=transcript)) + geom_line(stat="identity") + ggtitle("Genome Function at Fourth Evolution")

best5 <- read.table("best_replicated_gen_345.tsv", header=TRUE)
best5 <- filter(best5, species == "proteinX" | species == "proteinY" | species == "proteinZ")
best5_plot <-  ggplot(best5, aes(fill=species, color=species, x=time, y=transcript)) + geom_line(stat="identity") + ggtitle("Genome Function at Fifth Evolution")

best6 <- read.table("best_replicated_gen_348.tsv", header=TRUE)
best6 <- filter(best6, species == "proteinX" | species == "proteinY" | species == "proteinZ")
best6_plot <-  ggplot(best6, aes(fill=species, color=species, x=time, y=transcript)) + geom_line(stat="identity") + ggtitle("Genome Function at Sixth Evolution")

best7 <- read.table("best_replicated_gen_2590.tsv", header=TRUE)
best7 <- filter(best7, species == "proteinX" | species == "proteinY" | species == "proteinZ")
best7_plot <-  ggplot(best7, aes(fill=species, color=species, x=time, y=transcript)) + geom_line(stat="identity") + ggtitle("Genome Function at Seventh Evolution")

best8 <- read.table("best_replicated_gen_13635.tsv", header=TRUE)
best8 <- filter(best8, species == "proteinX" | species == "proteinY" | species == "proteinZ")
best8_plot <-  ggplot(best8, aes(fill=species, color=species, x=time, y=transcript)) + geom_line(stat="identity") + ggtitle("Genome Function at Eigth Evolution")

evo_plot <- grid.arrange(best1_plot, best2_plot, best3_plot, best4_plot, best5_plot, best6_plot, best7_plot, best8_plot, nrow=4, ncol=4)
ggsave("/home/wilkelab/pinetree-toys/three_genes_evolution/evolution_graphs.pdf", evo_plot, width=30, height=30, limitsize=FALSE)
