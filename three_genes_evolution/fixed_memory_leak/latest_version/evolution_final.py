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
ss_old = 600000
all_sos_list = [600000]
sos_iter_list = []
accepted = 0

#Setting arbitrary old fitness value
#eps = np.random.normal(mu, sigma)
#ss_new = ss_old * (1.0 + eps)

#Start of evolution program
while i < 1001:

    with open('new_gene.yml') as f:
        genome_tracker_old = yaml.safe_load(f)
    with open('new_gene.yml') as f:
        genome_tracker_new = yaml.safe_load(f)
    possibilities = [mutation_choices.modify_promoter(genome_tracker_new), mutation_choices.modify_rnase(genome_tracker_new), mutation_choices.modify_terminator(genome_tracker_new)]
    random.choice(possibilities)

    ss_new = mutation_test.test_mutation(df, genome_tracker_new, i)
    accept_prob = fitness_score.calc_fitness(ss_new, ss_old)
    if accept_prob > random.random():
        genome_tracker_old = genome_tracker_new
        save.save_file(pd.read_csv("three_genes_replicated.tsv", header=0, sep='\t'), "three_genes_replicated_{}.tsv".format(i))
        ss_old = ss_new
        all_sos_list.append(ss_new)
        sos_iter_list.append(i)
        accepted+=1
    else:
        pass

    i+=1
    print("i =", i)
f.close()

'''from mutation_accepted import gen, sos_list, all_sos_list, sos_iter_list, accepted
#Exported tsv files
all_sos_dataframe = pd.DataFrame(all_sos_list, columns=["Sum_of_Squares"])
export_csv = all_sos_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/fixed_memory_leak/latest_version/all_sos_data.tsv", index=False)'''
all_sos_list.remove(600000)
sos_dataframe = pd.DataFrame(data=all_sos_list, columns=["Sum_of_Squares"])
export_csv = sos_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/fixed_memory_leak/latest_version/all_sos_data.tsv", index=False)
sos_iter_dataframe = pd.DataFrame(data=sos_iter_list, columns=['Iteration'])
export_csv = sos_iter_dataframe.to_csv('~/pinetree-toys/three_genes_evolution/fixed_memory_leak/latest_version/sos_iter_data.tsv', index=False)
print("Accepted mutations =", accepted)
#print("Generations =", gen)
