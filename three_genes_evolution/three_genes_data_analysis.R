library(cowplot)
library(dplyr)

getwd()
setwd("/home/wilkelab/Projects/phage_simulation-master/data/simulation/three_genes")

data <- read.table("three_genes_rnase.tsv", header = TRUE)
data <- filter(data, species == c("proteinX", "proteinY", "proteinZ"), time > 99, time < 100)

box_plot1 <- ggplot(data, aes(fill=species, x=species, y=transcript)) + geom_bar(stat="identity")
box_plot1

data2 <- read.table("three_genes_rnase.tsv", header = TRUE)
data2 <- filter(data2, species == c("proteinX", "proteinY", "proteinZ"), time > 239, time < 240)

box_plot2 <- ggplot(data2, aes(fill=species, x=species, y=transcript)) + geom_bar(stat="identity")
box_plot2

data3 <- read.table("three_genes_rnase.tsv", header = TRUE)
data3 <- filter(data3, species == c("proteinX", "proteinY", "proteinZ"))

line_plot <- ggplot(data3, aes(fill=species, color=species, x=time, y=transcript)) + geom_point(stat="identity") + geom_line(stat="identity")
line_plot

setwd("/home/wilkelab/pinetree-toys/three_genes_evolution")
data4 <- read.table("all_sos_data.tsv", header=TRUE)
data4$index <- seq.int(nrow(data4))

sos_line_plot <- ggplot(data4, aes(x=index, y=Sum_of_Squares)) + geom_point(stat="identity") + geom_line(stat="identity")
sos_line_plot
