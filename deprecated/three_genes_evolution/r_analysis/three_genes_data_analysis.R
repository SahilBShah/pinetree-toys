library(cowplot)
library(dplyr)

getwd()
setwd("/home/sahil/pinetree-toys/three_genes_evolution/three_genes_target")

data <- read.table("three_genes_test_file1.tsv", header = TRUE)
data <- filter(data, species == "proteinX" | species == "proteinY" | species == "proteinZ", time > 98, time < 100)

box_plot1 <- ggplot(data, aes(fill=species, x=species, y=transcript)) + geom_bar(stat="identity") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
box_plot1

data2 <- read.table("three_genes_test_file1.tsv", header = TRUE)
data2 <- filter(data2, species == "proteinX" | species == "proteinY" | species == "proteinZ", time > 239, time < 240)

box_plot2 <- ggplot(data2, aes(fill=species, x=species, y=transcript)) + geom_bar(stat="identity") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
box_plot2

data3 <- read.table("three_genes_test_file.tsv", header = TRUE)
data3 <- filter(data3, species == "proteinX" | species == "proteinY" | species == "proteinZ")

line_plot <- ggplot(data3, aes(fill=species, color=species, x=time, y=transcript)) + geom_line(stat="identity") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
line_plot

data3 <- read.table("three_genes_replicated.tsv", header = TRUE)
data3 <- filter(data3, species == "proteinX" | species == "proteinY" | species == "proteinZ")

line_plot2 <- ggplot(data3, aes(fill=species, color=species, x=time, y=transcript)) + geom_point(stat="identity") + geom_line(stat="identity")
line_plot2

setwd("/home/wilkelab/pinetree-toys/three_genes_evolution")
data4 <- read.table("sos_data.tsv", header=TRUE)
data5 <- read.table("iter_data.tsv", header=TRUE)
data4 <- cbind(data4, data5)

sos_line_plot <- ggplot(data4, aes(x=Iteration, y=Sum_of_Squares)) + geom_point(stat="identity") + geom_line(stat="identity")
sos_line_plot
