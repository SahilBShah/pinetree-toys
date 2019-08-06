import pandas as pd
import string
import file_setup
import sum_of_squares
import save
import csv
import os_genome
import os

global gen
global sos_list
global all_sos_list
global sos_iter_list
gen = 0
sos_list = [60000]
all_sos_list = []
sos_iter_list = []

def accept_mutation(df, genome_tracker, i):
    """
    If the mutation is accepted based of the new fitness value or a calculated probability then the new genome will be simulated and saved if the sum of squares decreased.
    """

    global gen

    genome_tracker['output_file_name'] = "three_genes_replicated.tsv"
    genome_output_file_name = "genome_tracker_" + str(i) +".yml"

    #Accepting mutation
    os_genome.pt_call()
    #Taking in new file and removing unnecessary rows and columns
    nf = pd.read_csv("three_genes_replicated.tsv", header=0, sep='\t')
    saveable_nf = nf = pd.read_csv("three_genes_replicated.tsv", header=0, sep='\t')
    nf = file_setup.rearrange_file(nf)
    sos = sum_of_squares.calc_sum_of_squares(df, nf)
    all_sos_list.append(sos)
    #Accepts mutation only if sum of squares value decreases
    if sos <= sos_list[-1]:
        genome_tracker['output_file_name'] = "three_genes_replicated_" + str(i) + ".tsv"
        save.save_file(saveable_nf, output_file_name)
        save.save_file(genome_tracker, genome_output_file_name)
        sos_list.append(sos)
        sos_iter_list.append(i)
        gen+=1
