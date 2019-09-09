import pandas as pd
import string
import file_setup
import sum_of_squares
import save
import csv
from os_genome import call_pt
import yaml

global gen
global accepted
global sos_list
global all_sos_list
global sos_iter_list
gen = 0
accepted = 0
sos_list = [60000]
all_sos_list = []
sos_iter_list = []

def accept_mutation(df, genome_tracker, i):
    """
    If the mutation is accepted based of the new fitness value or a calculated probability then the new genome will be simulated and saved if the sum of squares decreased.
    """

    global gen
    global accepted

    genome_tracker['output_file_name'] = "three_genes_replicated.tsv"
    with open('new_gene.yml', 'w') as f:
        yaml.dump(genome_tracker, f)
    genome_output_file_name = "genome_tracker_" + str(i) +".yml"

    #Accepting mutation
    call_pt.pt_call()
    #Taking in new file and removing unnecessary rows and columns
    nf = pd.read_csv("three_genes_replicated.tsv", header=0, sep='\t')
    saveable_nf = nf = pd.read_csv("three_genes_replicated.tsv", header=0, sep='\t')
    nf = file_setup.rearrange_file(nf)
    sos = sum_of_squares.calc_sum_of_squares(df, nf)
    all_sos_list.append(sos)
    if i < 10001 and i > 9980:
        genome_tracker['output_file_name'] = "three_genes_replicated_" + str(i) + ".tsv"
        with open('new_gene.yml', 'w') as f:
            yaml.dump(genome_tracker, f)
        save.save_file(saveable_nf, genome_tracker['output_file_name'])
        with open(genome_output_file_name, 'w') as f:
            yaml.dump(genome_tracker, f)
    #Accepts mutation only if sum of squares value decreases
    if sos <= sos_list[-1]:
        genome_tracker['output_file_name'] = "three_genes_replicated_" + str(i) + ".tsv"
        with open('new_gene.yml', 'w') as f:
            yaml.dump(genome_tracker, f)
        save.save_file(saveable_nf, genome_tracker['output_file_name'])
        with open(genome_output_file_name, 'w') as f:
            yaml.dump(genome_tracker, f)
        sos_list.append(sos)
        sos_iter_list.append(i)
        gen+=1
    accepted+=1

    f.close()
