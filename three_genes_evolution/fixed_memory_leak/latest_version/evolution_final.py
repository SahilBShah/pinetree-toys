import evolution_setup
import file_setup
import save
import mutation_test
import mutation_accepted
import fitness_score
import sum_of_squares
import mutation_choices
import random
import numpy as np
import pandas as pd
import yaml

i = 1

input_files_list = evolution_setup.genome_setup()
df = input_files_list[0]
genome_tracker = input_files_list[1]
ss_old = 100000
all_sos_list = [100000]
sos_iter_list = []
accepted = 0


#Start of evolution program
while i < 101:

    with open('new_gene.yml') as f:
        genome_tracker_old = yaml.safe_load(f)
    with open('new_gene.yml') as f:
        genome_tracker_new = yaml.safe_load(f)
    possibilities = [mutation_choices.modify_promoter(genome_tracker_new), mutation_choices.modify_rnase(genome_tracker_new), mutation_choices.modify_terminator(genome_tracker_new)]
    random.choice(possibilities)

    ss_new = mutation_test.test_mutation(df, genome_tracker_new, i)
    accept_prob = fitness_score.calc_fitness(ss_new, ss_old)
    print("ss_new =", ss_new)
    if accept_prob > random.random():
        print("accepted")
        genome_tracker_old = genome_tracker_new
        #save.save_file(pd.read_csv("three_genes_replicated.tsv", header=0, sep='\t'), "three_genes_replicated_{}.tsv".format(i))
        ss_old = ss_new
        all_sos_list.append(ss_new)
        sos_iter_list.append(i)
        accepted+=1


    i+=1
    print("i =", i)


all_sos_list.remove(100000)
sos_dataframe = pd.DataFrame(data=all_sos_list, columns=["Sum_of_Squares"])
export_csv = sos_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/fixed_memory_leak/latest_version/all_sos_data.tsv", index=False)
sos_iter_dataframe = pd.DataFrame(data=sos_iter_list, columns=['Iteration'])
export_csv = sos_iter_dataframe.to_csv('~/pinetree-toys/three_genes_evolution/fixed_memory_leak/latest_version/sos_iter_data.tsv', index=False)
print("Accepted mutations =", accepted)
#print("Generations =", gen)
