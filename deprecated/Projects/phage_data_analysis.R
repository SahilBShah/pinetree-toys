library(cowplot)
library(dplyr)
library(readr)

setwd("/home/wilkelab/pinetree-toys/Projects/phage_simulation-master/data/simulation/phage/")
phage_data <- read_tsv('phage_wildtype.tsv')
phage_data <- dplyr::filter(phage_data, grepl('gene', species))
phage_data <- dplyr::filter(phage_data, !grepl('rbs', species))
phage_bar_plot <- ggplot(phage_data, aes(fill=species, x=species, y=transcript)) + geom_bar(stat="identity")
phage_bar_plot
