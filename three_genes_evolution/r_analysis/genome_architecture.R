library(cowplot)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(gtools)

setwd("/home/wilkelab/pinetree-toys/three_genes_evolution/genome_coordinates/best_genome")

gene1 <- read.table("gene1_data.tsv", header=TRUE)
gene2 <- read.table("gene2_data.tsv", header=TRUE)
gene3 <- read.table("gene3_data.tsv", header=TRUE)
promoter1 <- read.table("promoter1_data.tsv", header=TRUE)
promoter2 <- read.table("promoter2_data.tsv", header=TRUE)
terminator1 <- read.table("term1_data.tsv", header=TRUE)
terminator2 <- read.table("term2_data.tsv", header=TRUE)
terminator3 <- read.table("term3_data.tsv", header=TRUE)
rnase1 <- read.table("rnase1_data.tsv", header=TRUE)
rnase2 <- read.table("rnase2_data.tsv", header=TRUE)
rnase3 <- read.table("rnase3_data.tsv", header=TRUE)

y1 <- 1
y2 <- 2

if(promoter1[1, 1] > 0){
  draw_prom1 <- geom_segment(mapping=aes(x=20, y=2, xend=20, yend=2.5), arrow = arrow(length = unit(0.5, "cm")))
} else {
    draw_prom1 <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
}
if(promoter2[1, 1] > 0){
  draw_prom2 <- geom_segment(mapping=aes(x=40, y=2, xend=40, yend=2.5), arrow = arrow(length = unit(0.5, "cm")))
} else {
    draw_prom2 <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
}
if((terminator1[1, 1] >= 11) & (terminator1[2, 1] <= 26)){
  draw_term1 <- geom_segment(mapping=aes(x=10, y=2, xend=10, yend=2.25)) + geom_segment(mapping=aes(x=7, y=2.25, xend=13, yend=2.25))
} else {
    draw_term1 <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
}
if(terminator2[1, 1] > 0){
  draw_term2 <- geom_segment(mapping=aes(x=30, y=2, xend=30, yend=2.25)) + geom_segment(mapping=aes(x=27, y=2.25, xend=33, yend=2.25))
} else {
    draw_term2 <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
}
if(terminator3[1, 1] > 0){
  draw_term3 <- geom_segment(mapping=aes(x=50, y=2, xend=50, yend=2.25)) + geom_segment(mapping=aes(x=47, y=2.25, xend=53, yend=2.25))
} else {
    draw_term3 <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
}
if((rnase2[1, 1] >= 144) & (rnase2[2, 1] <= 159)){
  draw1_rnase2c <- geom_segment(mapping=aes(x=16, y=1.7, xend=20, yend=1.3))
  draw2_rnase2c <- geom_segment(mapping=aes(x=20, y=1.7, xend=16, yend=1.3))
} else {
    draw1_rnase2c <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
    draw2_rnase2c <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
}

genome_list <- list(draw_prom1, draw_prom2, draw_term1, draw_term2, draw_term3, draw1_rnase2c, draw2_rnase2c)

ggplot() + geom_rect(data=gene1, mapping=aes(xmin=1, xmax=10, ymin=y1, ymax=y2), fill="blue", alpha=0.5) + geom_rect(data=gene2, mapping=aes(xmin=20, xmax=30, ymin=y1, ymax=y2), fill="red", alpha=0.5) + geom_rect(data=gene3, mapping=aes(xmin=40, xmax=50, ymin=y1, ymax=y2), fill="yellow", alpha=0.5) + geom_segment(mapping=aes(x=1, y=2, xend=1, yend=2.5), col="black") + geom_segment(mapping=aes(x=30, y=1.5, xend=40, yend=1.5)) + geom_segment(mapping=aes(x=10, y=1.5, xend=20, yend=1.5), color="black") + geom_segment(mapping=aes(x=1, y=2.5, xend=3, yend=2.5), color="black", arrow = arrow(length = unit(0.5, "cm"))) + genome_list[1] + genome_list[2] + genome_list[3] + genome_list[4] + genome_list[5] + genome_list[6] + genome_list[7]

