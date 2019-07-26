import pandas as pd
import string
import file_setup
import sum_of_squares
import save
import genome_simulator
import csv

global gen
global sos_list
global all_sos_list
gen = 0
sos_list = [60000]
all_sos_list = []

def accept_mutation(df, genome_tracker, i):

    global gen

    output_file_name = "three_genes_replicated.tsv"
    genome_output_file_name = "genome_tracker_" + str(i) +".tsv"

    #Accepting mutation
    genome_simulator.simulate_genome(genome_tracker, output_file_name)
    #Taking in new file and removing unnecessary rows and columns
    nf = pd.read_csv("three_genes_replicated.tsv", header=0, sep='\t')
    nf = file_setup.rearrange_file(nf)
    sos = sum_of_squares.calc_sum_of_squares(df, nf)
    all_sos_list.append(sos)
    #Accepts mutation only if sum of squares value decreases
    if sos <= sos_list[-1]:
        output_file_name = "three_genes_replicated_" + str(i) + ".tsv"
        save.save_file(nf, output_file_name)
        save.save_file(genome_tracker, genome_output_file_name)
        sos_list.append(sos)
        gen+=1
