import pinetree as pt
import pandas as pd
import string
import random
import math
import itertools

#Add a genome_tracker return to baseline function

global Ne
global gen
global mu
global sigma
sos_list = [60000]
all_sos_list = []
f_old, genome_tracker['value'][12] = 1.0
mu = 0.0
sigma = 1.0
Ne = 10
i = 0
gen = 0

def main():

    df_name = input("Enter tsv file name: ")
    df = pd.read_csv(df_name, header=0, sep='\t')
    df = edit_target_file(df, df_name)
    genome_tracker = pd.read_csv("gene_tracker.tsv", header=0, sep='\t')

    #Setting random strengths for promoters and terminators
    genome_tracker['new_strength'][3], genome_tracker['previous_strength'][3] = random.randint(0, 3e10)
    genome_tracker['new_strength'][4], genome_tracker['previous_strength'][4] = random.randint(0, 3e10)
    genome_tracker['new_strength'][5], genome_tracker['previous_strength'][5] = random.randint(0, 3e10)
    genome_tracker['new_strength'][9], genome_tracker['previous_strength'][9] = random.uniform(0.0, 1.0)
    genome_tracker['new_strength'][10], genome_tracker['previous_strength'][10] = random.uniform(0.0, 1.0)
    genome_tracker['new_strength'][11], genome_tracker['previous_strength'][11] = random.uniform(0.0, 1.0)

    possibilities = ["alter polymerase 1 strength", "add promoter", "remove promoter", "add rnase",
                      "remove rnase", "add terminator", "remove terminator"]

    #Start of evolution program
    while i < 10000:

        mutation = random.choice(possibilities)

        if mutation == "alter polymerase 1 strength":
            promoter = "promoter"
            new_pol_strength = alter_poly_strength(genome_tracker, promoter)

        if mutation == "add promoter":
            new_promoter = add_promoter(genome_tracker)

        if mutation == "remove promoter":
            old_promoter = remove_promoter(genome_tracker)

        if mutation == "add rnase":
            new_rnase = add_rnase(genome_tracker)

        if mutation == "remove rnase":
            old_rnase = remove_rnase(genome_tracker)

        if mutation == "add terminator":
            new_terminator = add_terminator(genome_tracker)

        if mutation == "remove terminator":
            old_terminator = remove_terminator(genome_tracker)

        eps = np.random.normal(mu, sigma)
        f_new, genome_tracker['value'][13] = f_old * (1.0 + eps)
        test = accept_mutation(genome_tracker, i)

        if i == 0:
            output_file = "gen_0_data.tsv"
            three_genome.recreated_genome(genome_tracker, output_file)
            save_file(genome_tracker)

        if i == 4999:
            output_file = "gen_5000_data.tsv"
            three_genome.recreated_genome(genome_tracker, output_file)
            save_file(genome_tracker)

        if i == 9999:
            output_file = "gen_10000_data.tsv"
            three_genome.recreated_genome(genome_tracker, output_file)
            save_file(genome_tracker)

        i+=1
        print("i =", i)

    #Exported tsv files
    all_sos_dataframe = pd.DataFrame(all_sos_list, columns=["Sum_of_Squares"])
    export_csv = all_sos_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/all_sos_data.tsv", index=False)
    sos_list.remove(60000)
    sos_dataframe = pd.DataFrame(data=sos_list, columns=["Sum_of_Squares"])
    export_csv = sos_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/sos_data.tsv", index=False)

#Calculates the fitness of the new mutation
def calc_fitness(variant_fit, orig_fit, Ne):

    xi = orig_fit
    xj = variant_fit
    N = Ne

    if xj == xi:
        return(1.0 / float(N))
    if xj > xi:
        return(1.0)
    else:
        variant_fit = (xj / xi) ** (2 * N - 2)
        return variant_fit

    try:
        resolved = ((1-pow((xi/xj), 2)) / (1-pow((xi/xj), (2 * float(N)))))
    except OverflowError as e:
        resolved = 0.0
    return (resolved)

#Minimizes distance between observed values and regression line from given data
def calc_sum_of_squares(target_file, new_file):

    df = new_file
    df1 = target_file
    df_diff = df - df1
    df_squared = df_diff ** 2
    sum_df_squared = df_squared.sum()
    sos = sum_df_squared[0] + sum_df_squared[1] + sum_df_squared[2]
    return sos

#Removes unnecessary rows and columns in produced file
def edit_new_file(new_file):

    new_file = new_file[new_file['species'].isin(['proteinX', 'proteinY', 'proteinZ'])]
    new_file = new_file[['time', 'species', 'transcript']]
    new_file['time'] = new_file['time'].round().astype(int)
    new_file = new_file.pivot(index='time', columns='species', values='transcript')
    new_file = new_file.fillna(0.0)
    return new_file

#Removes unnecessary rows and columns in target file
def edit_target_file(target_file, name_of_file):

    target_file = target_file[target_file['species'].isin(['proteinX', 'proteinY', 'proteinZ'])]
    target_file = target_file[['time', 'species', 'transcript']]
    target_file = target_file.pivot(index='time', columns='species', values='transcript')
    target_file = target_file.fillna(0.0)
    return target_file

def accept_mutation(df, genome_tracker, i):

    global gen
    output_file = "three_genes_replicated.tsv" + str(i)

    if f_new > f_old:
        #Accepting mutation
        three_genome.recreated_genome(genome_tracker, output_file)
        #Taking in new file and removing unnecessary rows and columns
        nf = pd.read_csv("three_genes_replicated.tsv", header=0, sep='\t')
        nf = edit_new_file(nf)
        sos = sum_of_squares(df, nf)
        all_sos_list.append(sos)
        #Accepts mutation only if sum of squares value decreases
        if sos <= sos_list[-1]:
            save_file(nf)
            genome_ouput_file = "genome_tracker" + str(i)
            save_file(genome_coordinates, genome_output_file)
            sos_list.append(sos)
