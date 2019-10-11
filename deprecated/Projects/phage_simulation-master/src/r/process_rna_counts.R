###
# This script takes raw TSV tables of RNA-seq counts, normalizes them, maps IDs, and combines them into a single
# CSV containing all gene count information. 
###

library(tidyverse)

# Read in ID map
labels <- read_csv("../../output/id_map.csv") %>% 
  mutate(accession = trimws(accession)) # Trim an whitespace that causes problems with joins

import_counts <- function(file_name, replicate, time_point) {
  read_tsv(file_name, col_names=FALSE) %>%
    filter(X1 == "V01146.1", X3 == "CDS") %>%
    select("start" = X4, "end" = X5, "gene" = X9, "count" = X10) %>%
    mutate(gene_number = str_match(gene, "Note=(possible )?gene ([0-9\\.\\-AB]*)")[,3],
           lengthkb = (end - start)/1000,
           rpk = count/lengthkb,
           scale_factor = sum(rpk)/1000000,
           tpm = rpk/scale_factor) %>%  # Compute TPM
    left_join(labels, by = c("gene_number" = "gene_number"))
} 

###
# 37C time-course first
###
file_inputs <- bind_rows(
  list(file = "../../data/rna_seq/37C/mapped_to_T7_ecoli/PA36-wild10A-46-9min-rep1.counts.tsv", b_rep = 1, time = 9),
  list(file = "../../data/rna_seq/37C/mapped_to_T7_ecoli/PA45-wild10A-46-9min-rep2.counts.tsv", b_rep = 2, time = 9),
  list(file = "../../data/rna_seq/37C/mapped_to_T7_ecoli/PA54-wild10A-46-9min-rep3.counts.tsv", b_rep = 3, time = 9),
  list(file = "../../data/rna_seq/37C/mapped_to_T7_ecoli/PA35-wild10A-46-5min-rep1.counts.tsv", b_rep = 1, time = 5),
  list(file = "../../data/rna_seq/37C/mapped_to_T7_ecoli/PA44-wild10A-46-5min-rep2.counts.tsv", b_rep = 2, time = 5),
  list(file = "../../data/rna_seq/37C/mapped_to_T7_ecoli/PA53-wild10A-46-5min-rep3.counts.tsv", b_rep = 3, time = 5),
  list(file = "../../data/rna_seq/37C/mapped_to_T7_ecoli/PA34-wild10A-46-1min-rep1.counts.tsv", b_rep = 1, time = 1),
  list(file = "../../data/rna_seq/37C/mapped_to_T7_ecoli/PA43-wild10A-46-1min-rep2.counts.tsv", b_rep = 2, time = 1),
  list(file = "../../data/rna_seq/37C/mapped_to_T7_ecoli/PA52-wild10A-46-1min-rep3.counts.tsv", b_rep = 3, time = 1)
)

counts <- file_inputs %>% 
  mutate(read_counts = map(file, ~import_counts(.))) %>% # Run our import_counts function over all files
  unnest() %>%
  select(-file)

write_csv(counts, "../../output/transcripts_37C.csv")

###
# 25C time-course
###

file_inputs_25 <- bind_rows(
  list(file = "../../data/rna_seq/25C/mapped_to_T7_ecoli/PA196_T7W-R1-05m_S27_L004_R1_001_sort_counts.tsv", b_rep = 1, time = 5),
  list(file = "../../data/rna_seq/25C/mapped_to_T7_ecoli/PA206_T7W-R2-05m_S32_L004_R1_001_sort_counts.tsv", b_rep = 2, time = 5),
  list(file = "../../data/rna_seq/25C/mapped_to_T7_ecoli/PA216_T7W-R3-05m_S37_L004_R1_001_sort_counts.tsv", b_rep = 3, time = 5),
  list(file = "../../data/rna_seq/25C/mapped_to_T7_ecoli/PA198_T7W-R1-10m_S28_L004_R1_001_sort_counts.tsv", b_rep = 1, time = 10),
  list(file = "../../data/rna_seq/25C/mapped_to_T7_ecoli/PA208_T7W-R2-10m_S33_L004_R1_001_sort_counts.tsv", b_rep = 2, time = 10),
  list(file = "../../data/rna_seq/25C/mapped_to_T7_ecoli/PA218_T7W-R3-10m_S38_L004_R1_001_sort_counts.tsv", b_rep = 3, time = 10),
  list(file = "../../data/rna_seq/25C/mapped_to_T7_ecoli/PA200_T7W-R1-15m_S29_L004_R1_001_sort_counts.tsv", b_rep = 1, time = 15),
  list(file = "../../data/rna_seq/25C/mapped_to_T7_ecoli/PA210_T7W-R2-15m_S34_L004_R1_001_sort_counts.tsv", b_rep = 2, time = 15),
  list(file = "../../data/rna_seq/25C/mapped_to_T7_ecoli/PA220_T7W-R3-15m_S39_L004_R1_001_sort_counts.tsv", b_rep = 3, time = 15),
  list(file = "../../data/rna_seq/25C/mapped_to_T7_ecoli/PA202_T7W-R1-20m_S30_L004_R1_001_sort_counts.tsv", b_rep = 1, time = 20),
  list(file = "../../data/rna_seq/25C/mapped_to_T7_ecoli/PA212_T7W-R2-20m_S35_L004_R1_001_sort_counts.tsv", b_rep = 2, time = 20),
  list(file = "../../data/rna_seq/25C/mapped_to_T7_ecoli/PA222_T7W-R3-20m_S40_L004_R1_001_sort_counts.tsv", b_rep = 3, time = 20),
  list(file = "../../data/rna_seq/25C/mapped_to_T7_ecoli/PA204_T7W-R1-25m_S31_L004_R1_001_sort_counts.tsv", b_rep = 1, time = 25),
  list(file = "../../data/rna_seq/25C/mapped_to_T7_ecoli/PA214_T7W-R2-25m_S36_L004_R1_001_sort_counts.tsv", b_rep = 2, time = 25),
  list(file = "../../data/rna_seq/25C/mapped_to_T7_ecoli/PA224_T7W-R3-25m_S41_L004_R1_001_sort_counts.tsv", b_rep = 3, time = 25)
)

counts_25 <- file_inputs_25 %>% 
  mutate(read_counts = map(file, ~import_counts(.))) %>% # Run our import_counts function over all files
  unnest() %>%
  select(-file)

write_csv(counts_25, "../../output/transcripts_25C.csv")
