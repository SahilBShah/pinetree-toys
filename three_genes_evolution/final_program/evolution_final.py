import evolution_setup
import file_setup
import mutation_test
import fitness_score
import sum_of_squares
import mutation_choices
import random
import numpy as np
import pandas as pd
import yaml
import initialize_yaml

#General setup
ss_old = 10000000
all_sos_list = [ss_old]
sos_iter_list = []
is_accepted = []
accepted = 0
i = 1


#Target
input_df = input("Enter tsv file name: ")
df = pd.read_csv(input_df, header=0, sep='\t')
df = file_setup.rearrange_file(df)
#df = evolution_setup.genome_setup()


#ddhslgkash Opens yaml files containing genome coordinates
initialize_yaml.create_yaml()
with open('new_gene.yml') as f:
    genome_tracker_old = yaml.safe_load(f)
with open('new_gene.yml') as f:
    genome_tracker_new = yaml.safe_load(f)

#Start of evolution program
while i < 5001:

    #Mutation is chosen and performed on best genome
    possibilities = [mutation_choices.modify_promoter, mutation_choices.modify_rnase, mutation_choices.modify_terminator]
    random.choice(possibilities)(genome_tracker_new)
    with open('new_gene.yml') as f:
        genome_tracker_new = yaml.safe_load(f)

    #Sum of squares is calculated and the mutation is accepted or rejected based off of its calculated fitness value
    ss_new = mutation_test.test_mutation(df)
    accept_prob = fitness_score.calc_fitness(ss_new, ss_old)
    all_sos_list.append(ss_new)
    sos_iter_list.append(i)
    if accept_prob > random.random():
        #If accepted....
        genome_tracker_old = genome_tracker_new
        with open('gene_{}.yml'.format(i), 'w') as f:
            yaml.dump(genome_tracker_old, f)
        save_df = pd.read_csv("three_genes_replicated.tsv", header=0, sep='\t')
        save_df.to_csv("three_genes_replicated_{}.tsv".format(i), sep='\t', index=False)
        ss_old = ss_new
        is_accepted.append("yes")
        accepted+=1
    else:
        #a bit about me
        genome_tracker_new = genome_tracker_old
        is_accepted.append("no")


    i+=1
    print("i =", i)


all_sos_list = all_sos_list[1:]
sos_dataframe = pd.DataFrame(data=zip(sos_iter_list, all_sos_list, is_accepted), columns=["Iteration", "Sum_of_Squares", "Accepted"])
export_csv = sos_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/final_program/sos_data.tsv", index=False, sep='\t')
print("Accepted mutations =", accepted)
