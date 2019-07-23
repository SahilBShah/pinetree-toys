import pandas as pd
import csv
import file_setup
import random

def genome_setup():

  input_df = input("Enter tsv file name: ")
  df = pd.read_csv(input_df, header=0, sep='\t')
  df = file_setup.rearrange_file(df)
  genome_tracker = pd.read_csv("gene_tracker.tsv", header=0, sep='\t')
  genome_tracker.reset_index()
  genome_tracker.set_index('species', inplace=True)
  genome_tracker.loc[('f_old'), ('value')] = 1.0

  #Setting random strengths for promoters and terminators
  genome_tracker.loc[('promoter1'), ('previous_strength')] = genome_tracker.loc[('promoter1'), ('current_strength')]
  genome_tracker.loc[('promoter1'), ('current_strength')] = random.randint(0, 3e10)
  return [df, genome_tracker]
