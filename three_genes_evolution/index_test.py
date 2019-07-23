import pandas as pd
import csv

genome_tracker = pd.read_csv('gene_tracker.tsv', header=0, sep='\t')
genome_tracker.reset_index()
genome_tracker.set_index('species', inplace=True)
a = genome_tracker['start']['promoter1'] = 10
#print(a)
