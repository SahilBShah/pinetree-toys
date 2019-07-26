import evolution_setup
import file_setup
import save
import mutation_test
import mutation_accepted
import fitness_score
import sum_of_squares
import mutation_choices
import genome_simulator
import random
import numpy as np
import pandas as pd

mu = 0.0
sigma = 1.0
i = 1

input_files_list = evolution_setup.genome_setup()
df = input_files_list[0]
genome_tracker = input_files_list[1]

#Start of evolution program
while i < 10001:

    possibilities = [mutation_choices.modify_promoter(genome_tracker), mutation_choices.modify_rnase(genome_tracker), mutation_choices.modify_terminator(genome_tracker)]
    random.choice(possibilities)

    eps = np.random.normal(mu, sigma)
    genome_tracker.loc[('f_new'), ('value')] = genome_tracker.loc[('f_old'), ('value')] * (1.0 + eps)
    mutation_test.test_mutation(df, genome_tracker, i)

    i+=1
    print("i =", i)

from mutation_accepted import gen, sos_list, all_sos_list
#Exported tsv files
all_sos_dataframe = pd.DataFrame(all_sos_list, columns=["Sum_of_Squares"])
export_csv = all_sos_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/all_sos_data.tsv", index=False)
sos_list.remove(60000)
sos_dataframe = pd.DataFrame(data=sos_list, columns=["Sum_of_Squares"])
export_csv = sos_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/sos_data.tsv", index=False)
print("Generations =", gen)
