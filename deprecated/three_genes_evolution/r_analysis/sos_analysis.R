library(cowplot)
library(dplyr)

setwd("/home/wilkelab/pinetree-toys/three_genes_evolution")
data4 <- read.table("sos_data.tsv", header=TRUE)
data5 <- read.table("iter_data.tsv", header=TRUE)
data4 <- cbind(data4, data5)

sos_line_plot <- ggplot(data4, aes(x=Iteration, y=Sum_of_Squares)) + geom_point(stat="identity") + geom_line(stat="identity")
sos_line_plot