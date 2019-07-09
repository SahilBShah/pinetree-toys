library(cowplot)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(gtools)

#Final genome data
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

#Draw promoters
if((promoter1[1, 1] >= 122) & (promoter1[2, 1] <= 132)){
  draw_prom1_r1a <- geom_segment(mapping=aes(x=12, y=1.5, xend=12, yend=2))
  draw_prom1_r2a <- geom_segment(mapping=aes(x=12, y=2, xend=13, yend=2), arrow = arrow(length = unit(0.5, "cm")))
} else {
  draw_prom1_r1a <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
  draw_prom1_r2a <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
}
if((promoter1[1, 1] >= 133) & (promoter1[2, 1] <= 143)){
  draw_prom1_r1b <- geom_segment(mapping=aes(x=15, y=1.5, xend=15, yend=2))
  draw_prom1_r2b <- geom_segment(mapping=aes(x=15, y=2, xend=16, yend=2), arrow = arrow(length = unit(0.5, "cm")))
} else {
  draw_prom1_r1b <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
  draw_prom1_r2b <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
}
if((promoter2[1, 1] >= 281) & promoter2[2, 1] <= 291){
  draw_prom2_r1a <- geom_segment(mapping=aes(x=32, y=1.5, xend=32, yend=2))
  draw_prom2_r2a <- geom_segment(mapping=aes(x=32, y=2, xend=33, yend=2), arrow = arrow(length = unit(0.5, "cm")))
} else {
  draw_prom2_r1a <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
  draw_prom2_r2a <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
}
if((promoter2[1, 1] >= 292) & promoter2[2, 1] <= 302){
  draw_prom2_r1b <- geom_segment(mapping=aes(x=35, y=1.5, xend=35, yend=2))
  draw_prom2_r2b <- geom_segment(mapping=aes(x=35, y=2, xend=36, yend=2), arrow = arrow(length = unit(0.5, "cm")))
} else {
  draw_prom2_r1b <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
  draw_prom2_r2b <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
}

#Draw terminators
if((terminator1[1, 1] >= 122) & (terminator1[2, 1] <= 132)){
  draw_term1_r1a <- geom_segment(mapping=aes(x=12, y=1.5, xend=12, yend=1.75))
  draw_term1_r2a <- geom_segment(mapping=aes(x=11, y=1.75, xend=13, yend=1.75))
} else {
  draw_term1_r1a <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
  draw_term1_r2a <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
}
if((terminator1[1, 1] >= 133) & (terminator1[2, 1] <= 143)){
  draw_term1_r1b <- geom_segment(mapping=aes(x=15, y=1.5, xend=15, yend=1.75))
  draw_term1_r2b <- geom_segment(mapping=aes(x=14, y=1.75, xend=16, yend=1.75))
} else {
  draw_term1_r1b <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
  draw_term1_r2b <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
}
if((terminator1[1, 1] >= 144) & (terminator1[2, 1] <= 159)){
  draw_term1_r1c <- geom_segment(mapping=aes(x=18, y=1.5, xend=18, yend=1.75))
  draw_term1_r2c <- geom_segment(mapping=aes(x=17, y=1.75, xend=19, yend=1.75))
} else {
  draw_term1_r1c <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
  draw_term1_r2c <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
}
if((terminator2[1, 1] >= 281) & (terminator2[2, 1] <= 291)){
  draw_term2_r1a <- geom_segment(mapping=aes(x=32, y=1.5, xend=32, yend=1.75))
  draw_term2_r2a <- geom_segment(mapping=aes(x=31, y=1.75, xend=33, yend=1.75))
} else {
  draw_term2_r1a <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
  draw_term2_r2a <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
}
if((terminator2[1, 1] >= 292) & (terminator2[2, 1] <= 302)){
  draw_term2_r1b <- geom_segment(mapping=aes(x=35, y=1.5, xend=35, yend=1.75))
  draw_term2_r2b <- geom_segment(mapping=aes(x=34, y=1.75, xend=36, yend=1.75))
} else {
  draw_term2_r1b <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
  draw_term2_r2b <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
}
if((terminator2[1, 1] >= 303) & (terminator2[2, 1] <= 318)){
  draw_term2_r1c <- geom_segment(mapping=aes(x=38, y=1.5, xend=38, yend=1.75))
  draw_term2_r2c <- geom_segment(mapping=aes(x=37, y=1.75, xend=39, yend=1.75))
} else {
  draw_term2_r1c <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
  draw_term2_r2c <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
}
if(terminator3[1, 1] > 0){
  draw_term3a <- geom_segment(mapping=aes(x=50, y=2, xend=50, yend=2.25))
  draw_term3b <- geom_segment(mapping=aes(x=47, y=2.25, xend=53, yend=2.25))
} else {
    draw_term3a <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
    draw_term3b <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
}

#Draw rnases (first segment, region a within first rnase placement)
if((rnase1[1, 1] >= 11) & (rnase1[2, 1] <= 25)){
  draw_rnase1a <- geom_segment(mapping=(aes(x=1, y=2.1, xend=3, yend=2)))
  draw_rnase1b <- geom_segment(mapping=(aes(x=3, y=2.1, xend=1, yend=2)))
}else{
  draw_rnase1a <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
  draw_rnase1b <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
}
if((rnase2[1, 1] >= 122) & (rnase2[1, 1] <= 132)){
  draw_rnase2_r1a <- geom_segment(mapping=aes(x=11, y=1.7, xend=13, yend=1.3))
  draw_rnase2_r2a <- geom_segment(mapping=aes(x=13, y=1.7, xend=11, yend=1.3))
}else{
  draw_rnase2_r1a <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
  draw_rnase2_r2a <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
}
if((rnase2[1, 1] >= 133) & (rnase2[1, 1] <= 143)){
  draw_rnase2_r1b <- geom_segment(mapping=aes(x=14, y=1.7, xend=16, yend=1.3))
  draw_rnase2_r2b <- geom_segment(mapping=aes(x=16, y=1.7, xend=14, yend=1.3))
}else{
  draw_rnase2_r1b <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
  draw_rnase2_r2b <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
}
if((rnase2[1, 1] >= 144) & (rnase2[2, 1] <= 159)){
  draw_rnase2_r1c <- geom_segment(mapping=aes(x=17, y=1.7, xend=19, yend=1.3))
  draw_rnase2_r2c <- geom_segment(mapping=aes(x=19, y=1.7, xend=17, yend=1.3))
} else {
    draw_rnase2_r1c <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
    draw_rnase2_r2c <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
}
if((rnase3[1, 1] >= 281) & (rnase3[1, 1] <= 291)){
  draw_rnase3_r1a <- geom_segment(mapping=aes(x=31, y=1.7, xend=33, yend=1.3))
  draw_rnase3_r2a <- geom_segment(mapping=aes(x=33, y=1.7, xend=31, yend=1.3))
}else{
  draw_rnase3_r1a <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
  draw_rnase3_r2a <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
}
if((rnase3[1, 1] >= 292) & (rnase3[1, 1] <= 302)){
  draw_rnase3_r1b <- geom_segment(mapping=aes(x=34, y=1.7, xend=36, yend=1.3))
  draw_rnase3_r2b <- geom_segment(mapping=aes(x=36, y=1.7, xend=34, yend=1.3))
}else{
  draw_rnase3_r1b <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
  draw_rnase3_r2b <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
}
if((rnase3[1, 1] >= 303) & (rnase3[1, 1] <= 318)){
  draw_rnase3_r1c <- geom_segment(mapping=aes(x=37, y=1.7, xend=39, yend=1.3))
  draw_rnase3_r2c <- geom_segment(mapping=aes(x=39, y=1.7, xend=37, yend=1.3))
}else{
  draw_rnase3_r1c <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
  draw_rnase3_r2c <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
}

genome_list <- list(draw_prom1_r1a, draw_prom1_r2a, draw_prom1_r1b, draw_prom1_r2b, draw_prom2_r1a, draw_prom2_r2a, draw_prom2_r1b, draw_prom2_r2b, draw_term1_r1a, draw_term1_r2a, draw_term1_r1b, draw_term1_r2b, draw_term1_r1c, draw_term1_r2c, draw_term2_r1a, draw_term2_r2a, draw_term2_r1b, draw_term2_r2b, draw_term2_r1c, draw_term2_r2c, draw_term3a, draw_term3b, draw_rnase1a, draw_rnase1b, draw_rnase2_r1a, draw_rnase2_r2a, draw_rnase2_r1b, draw_rnase2_r2b, draw_rnase2_r1c, draw_rnase2_r2c, draw_rnase3_r1a, draw_rnase3_r2a, draw_rnase3_r1b, draw_rnase3_r2b, draw_rnase3_r1c, draw_rnase3_r2c)

final_genome_plot <- ggplot() + ggtitle("Final Genome Architecture") + geom_rect(data=gene1, mapping=aes(xmin=1, xmax=10, ymin=y1, ymax=y2), fill="blue", alpha=0.5) + geom_rect(data=gene2, mapping=aes(xmin=20, xmax=30, ymin=y1, ymax=y2), fill="red", alpha=0.5) + geom_rect(data=gene3, mapping=aes(xmin=40, xmax=50, ymin=y1, ymax=y2), fill="yellow", alpha=0.5) + geom_segment(mapping=aes(x=1, y=2, xend=1, yend=2.5), col="black") + geom_segment(mapping=aes(x=30, y=1.5, xend=40, yend=1.5)) + geom_segment(mapping=aes(x=10, y=1.5, xend=20, yend=1.5), color="black") + geom_segment(mapping=aes(x=1, y=2.5, xend=3, yend=2.5), color="black", arrow = arrow(length = unit(0.5, "cm"))) + draw_prom1_r1a + draw_prom1_r2a + draw_prom1_r1b + draw_prom1_r2b + draw_prom2_r1a + draw_prom2_r2a + draw_prom2_r1b + draw_prom2_r2b + draw_term1_r1a + draw_term1_r2a + draw_term1_r1b + draw_term1_r2b + draw_term1_r1c + draw_term1_r2c + draw_term2_r1a + draw_term2_r2a + draw_term2_r1b + draw_term2_r2b + draw_term2_r1c + draw_term2_r2c + draw_term3a + draw_term3b + draw_rnase1a + draw_rnase1b + draw_rnase2_r1a + draw_rnase2_r2a + draw_rnase2_r1b + draw_rnase2_r2b + draw_rnase2_r1c + draw_rnase2_r2c + draw_rnase3_r1a + draw_rnase3_r2a + draw_rnase3_r1b + draw_rnase3_r2b + draw_rnase3_r1c + draw_rnase3_r2c



#Target genome data
setwd("/home/wilkelab/pinetree-toys/three_genes_evolution/genome_coordinates/target_genome")

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

#Draw promoters
if((promoter1[1, 1] >= 122) & (promoter1[2, 1] <= 132)){
  draw_prom1_r1a <- geom_segment(mapping=aes(x=12, y=1.5, xend=12, yend=2))
  draw_prom1_r2a <- geom_segment(mapping=aes(x=12, y=2, xend=13, yend=2), arrow = arrow(length = unit(0.5, "cm")))
} else {
  draw_prom1_r1a <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
  draw_prom1_r2a <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
}
if((promoter1[1, 1] >= 133) & (promoter1[2, 1] <= 143)){
  draw_prom1_r1b <- geom_segment(mapping=aes(x=15, y=1.5, xend=15, yend=2))
  draw_prom1_r2b <- geom_segment(mapping=aes(x=15, y=2, xend=16, yend=2), arrow = arrow(length = unit(0.5, "cm")))
} else {
  draw_prom1_r1b <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
  draw_prom1_r2b <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
}
if((promoter2[1, 1] >= 281) & promoter2[2, 1] <= 291){
  draw_prom2_r1a <- geom_segment(mapping=aes(x=32, y=1.5, xend=32, yend=2))
  draw_prom2_r2a <- geom_segment(mapping=aes(x=32, y=2, xend=33, yend=2), arrow = arrow(length = unit(0.5, "cm")))
} else {
  draw_prom2_r1a <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
  draw_prom2_r2a <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
}
if((promoter2[1, 1] >= 292) & promoter2[2, 1] <= 302){
  draw_prom2_r1b <- geom_segment(mapping=aes(x=35, y=1.5, xend=35, yend=2))
  draw_prom2_r2b <- geom_segment(mapping=aes(x=35, y=2, xend=36, yend=2), arrow = arrow(length = unit(0.5, "cm")))
} else {
  draw_prom2_r1b <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
  draw_prom2_r2b <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
}

#Draw terminators
if((terminator1[1, 1] >= 122) & (terminator1[2, 1] <= 132)){
  draw_term1_r1a <- geom_segment(mapping=aes(x=12, y=1.5, xend=12, yend=1.75))
  draw_term1_r2a <- geom_segment(mapping=aes(x=11, y=1.75, xend=13, yend=1.75))
} else {
  draw_term1_r1a <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
  draw_term1_r2a <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
}
if((terminator1[1, 1] >= 133) & (terminator1[2, 1] <= 143)){
  draw_term1_r1b <- geom_segment(mapping=aes(x=15, y=1.5, xend=15, yend=1.75))
  draw_term1_r2b <- geom_segment(mapping=aes(x=14, y=1.75, xend=16, yend=1.75))
} else {
  draw_term1_r1b <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
  draw_term1_r2b <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
}
if((terminator1[1, 1] >= 144) & (terminator1[2, 1] <= 159)){
  draw_term1_r1c <- geom_segment(mapping=aes(x=18, y=1.5, xend=18, yend=1.75))
  draw_term1_r2c <- geom_segment(mapping=aes(x=17, y=1.75, xend=19, yend=1.75))
} else {
  draw_term1_r1c <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
  draw_term1_r2c <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
}
if((terminator2[1, 1] >= 281) & (terminator2[2, 1] <= 291)){
  draw_term2_r1a <- geom_segment(mapping=aes(x=32, y=1.5, xend=32, yend=1.75))
  draw_term2_r2a <- geom_segment(mapping=aes(x=31, y=1.75, xend=33, yend=1.75))
} else {
  draw_term2_r1a <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
  draw_term2_r2a <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
}
if((terminator2[1, 1] >= 292) & (terminator2[2, 1] <= 302)){
  draw_term2_r1b <- geom_segment(mapping=aes(x=35, y=1.5, xend=35, yend=1.75))
  draw_term2_r2b <- geom_segment(mapping=aes(x=34, y=1.75, xend=36, yend=1.75))
} else {
  draw_term2_r1b <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
  draw_term2_r2b <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
}
if((terminator2[1, 1] >= 303) & (terminator2[2, 1] <= 318)){
  draw_term2_r1c <- geom_segment(mapping=aes(x=38, y=1.5, xend=38, yend=1.75))
  draw_term2_r2c <- geom_segment(mapping=aes(x=37, y=1.75, xend=39, yend=1.75))
} else {
  draw_term2_r1c <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
  draw_term2_r2c <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
}
if(terminator3[1, 1] > 0){
  draw_term3a <- geom_segment(mapping=aes(x=50, y=2, xend=50, yend=2.25))
  draw_term3b <- geom_segment(mapping=aes(x=47, y=2.25, xend=53, yend=2.25))
} else {
  draw_term3a <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
  draw_term3b <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
}

#Draw rnases (first segment, region a within first rnase placement)
if((rnase1[1, 1] >= 11) & (rnase1[2, 1] <= 25)){
  draw_rnase1a <- geom_segment(mapping=(aes(x=1, y=2.1, xend=3, yend=2)))
  draw_rnase1b <- geom_segment(mapping=(aes(x=3, y=2.1, xend=1, yend=2)))
}else{
  draw_rnase1a <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
  draw_rnase1b <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
}
if((rnase2[1, 1] >= 122) & (rnase2[1, 1] <= 132)){
  draw_rnase2_r1a <- geom_segment(mapping=aes(x=11, y=1.7, xend=13, yend=1.3))
  draw_rnase2_r2a <- geom_segment(mapping=aes(x=13, y=1.7, xend=11, yend=1.3))
}else{
  draw_rnase2_r1a <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
  draw_rnase2_r2a <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
}
if((rnase2[1, 1] >= 133) & (rnase2[1, 1] <= 143)){
  draw_rnase2_r1b <- geom_segment(mapping=aes(x=14, y=1.7, xend=16, yend=1.3))
  draw_rnase2_r2b <- geom_segment(mapping=aes(x=16, y=1.7, xend=14, yend=1.3))
}else{
  draw_rnase2_r1b <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
  draw_rnase2_r2b <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
}
if((rnase2[1, 1] >= 144) & (rnase2[2, 1] <= 159)){
  draw_rnase2_r1c <- geom_segment(mapping=aes(x=17, y=1.7, xend=19, yend=1.3))
  draw_rnase2_r2c <- geom_segment(mapping=aes(x=19, y=1.7, xend=17, yend=1.3))
} else {
  draw_rnase2_r1c <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
  draw_rnase2_r2c <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
}
if((rnase3[1, 1] >= 281) & (rnase3[1, 1] <= 291)){
  draw_rnase3_r1a <- geom_segment(mapping=aes(x=31, y=1.7, xend=33, yend=1.3))
  draw_rnase3_r2a <- geom_segment(mapping=aes(x=33, y=1.7, xend=31, yend=1.3))
}else{
  draw_rnase3_r1a <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
  draw_rnase3_r2a <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
}
if((rnase3[1, 1] >= 292) & (rnase3[1, 1] <= 302)){
  draw_rnase3_r1b <- geom_segment(mapping=aes(x=34, y=1.7, xend=36, yend=1.3))
  draw_rnase3_r2b <- geom_segment(mapping=aes(x=36, y=1.7, xend=34, yend=1.3))
}else{
  draw_rnase3_r1b <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
  draw_rnase3_r2b <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
}
if((rnase3[1, 1] >= 303) & (rnase3[1, 1] <= 318)){
  draw_rnase3_r1c <- geom_segment(mapping=aes(x=37, y=1.7, xend=39, yend=1.3))
  draw_rnase3_r2c <- geom_segment(mapping=aes(x=39, y=1.7, xend=37, yend=1.3))
}else{
  draw_rnase3_r1c <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
  draw_rnase3_r2c <- geom_segment(mapping=aes(x=0, y=0, xend=0, yend=0))
}

original_genome_plot <- ggplot() + ggtitle("Original Genome Architecture") + geom_rect(data=gene1, mapping=aes(xmin=1, xmax=10, ymin=y1, ymax=y2), fill="blue", alpha=0.5) + geom_rect(data=gene2, mapping=aes(xmin=20, xmax=30, ymin=y1, ymax=y2), fill="red", alpha=0.5) + geom_rect(data=gene3, mapping=aes(xmin=40, xmax=50, ymin=y1, ymax=y2), fill="yellow", alpha=0.5) + geom_segment(mapping=aes(x=1, y=2, xend=1, yend=2.5), col="black") + geom_segment(mapping=aes(x=30, y=1.5, xend=40, yend=1.5)) + geom_segment(mapping=aes(x=10, y=1.5, xend=20, yend=1.5), color="black") + geom_segment(mapping=aes(x=1, y=2.5, xend=3, yend=2.5), color="black", arrow = arrow(length = unit(0.5, "cm"))) + draw_prom1_r1a + draw_prom1_r2a + draw_prom1_r1b + draw_prom1_r2b + draw_prom2_r1a + draw_prom2_r2a + draw_prom2_r1b + draw_prom2_r2b + draw_term1_r1a + draw_term1_r2a + draw_term1_r1b + draw_term1_r2b + draw_term1_r1c + draw_term1_r2c + draw_term2_r1a + draw_term2_r2a + draw_term2_r1b + draw_term2_r2b + draw_term2_r1c + draw_term2_r2c + draw_term3a + draw_term3b + draw_rnase1a + draw_rnase1b + draw_rnase2_r1a + draw_rnase2_r2a + draw_rnase2_r1b + draw_rnase2_r2b + draw_rnase2_r1c + draw_rnase2_r2c + draw_rnase3_r1a + draw_rnase3_r2a + draw_rnase3_r1b + draw_rnase3_r2b + draw_rnase3_r1c + draw_rnase3_r2c

genome_plot <- grid.arrange(original_genome_plot, final_genome_plot, nrow=2, ncol=1)

ggsave("/home/wilkelab/pinetree-toys/three_genes_evolution/genome_coordinates/best_genome/genome_plot.pdf", genome_plot, width=30, height=30, limitsize=FALSE)

