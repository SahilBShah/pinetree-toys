import pinetree as pt
import pandas as pd
import string
import random
import math
import itertools
import numpy as np

#Add a genome_tracker return to baseline function
def main():

    global Ne
    global gen
    global mu
    global sigma
    global all_sos_list
    global sos_list
    sos_list = [60000]
    all_sos_list = []
    mu = 0.0
    sigma = 1.0
    Ne = 10
    i = 0
    gen = 0

    df_name = input("Enter tsv file name: ")
    df = pd.read_csv(df_name, header=0, sep='\t')
    df = edit_target_file(df, df_name)
    genome_tracker = pd.read_csv("gene_tracker.tsv", header=0, sep='\t')
    genome_tracker['value'].values[12] = 1.0

    #Setting random strengths for promoters and terminators
    genome_tracker['new_strength'].values[3] = genome_tracker['previous_strength'].values[3] = random.randint(0, 3e10)

    possibilities = ["alter polymerase 1 strength", "add promoter", "remove promoter", "add rnase",
                      "remove rnase", "add terminator", "remove terminator"]

    #Start of evolution program
    while i < 1000:

        mutation = random.choice(possibilities)

        if mutation == "alter polymerase 1 strength":
            which_promoter = "promoter"
            alter_poly_strength(genome_tracker, which_promoter)

        if mutation == "add promoter":
            add_promoter(genome_tracker)

        if mutation == "remove promoter":
            remove_promoter(genome_tracker)

        if mutation == "add rnase":
            add_rnase(genome_tracker)

        if mutation == "remove rnase":
            remove_rnase(genome_tracker)

        if mutation == "add terminator":
            add_terminator(genome_tracker)

        if mutation == "remove terminator":
            remove_terminator(genome_tracker)

        eps = np.random.normal(mu, sigma)
        genome_tracker['value'].values[13] = genome_tracker['value'].values[12] * (1.0 + eps)
        accept_mutation(df, genome_tracker, i)

        if i == 0:
            output_file_name = "gen_0_data.tsv"
            genome_output_file_name = "genome_tracker_" + str(i) + ".tsv"
            simulate_genome(genome_tracker, output_file_name)
            save_file(genome_tracker, genome_output_file_name)

        if i == 4999:
            output_file_name = "gen_5000_data.tsv"
            genome_output_file_name = "genome_tracker_" + str(i) + ".tsv"
            simulate_genome(genome_tracker, output_file_name)
            save_file(genome_tracker, genome_output_file_name)

        if i == 9999:
            output_file_name = "gen_10000_data.tsv"
            genome_output_file_name = "genome_tracker_" + str(i) + ".tsv"
            simulate_genome(genome_tracker, output_file_name)
            save_file(genome_tracker, genome_output_file_name)

        i+=1
        print("i =", i)

    #Exported tsv files
    all_sos_dataframe = pd.DataFrame(all_sos_list, columns=["Sum_of_Squares"])
    export_csv = all_sos_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/all_sos_data.tsv", index=False)
    sos_list.remove(60000)
    sos_dataframe = pd.DataFrame(data=sos_list, columns=["Sum_of_Squares"])
    export_csv = sos_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/sos_data.tsv", index=False)
    print("Generations =", gen)

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
    output_file_name = "three_genes_replicated.tsv"
    genome_output_file_name = "genome_tracker_" + str(i) +".tsv"
    f_new = genome_tracker['value'].values[13]
    f_old = genome_tracker['value'].values[12]

    #f_new > f_old
    if f_new > f_old:
        #Accepting mutation
        simulate_genome(genome_tracker, output_file_name)
        #Taking in new file and removing unnecessary rows and columns
        nf = pd.read_csv("three_genes_replicated.tsv", header=0, sep='\t')
        nf = edit_new_file(nf)
        sos = calc_sum_of_squares(df, nf)
        all_sos_list.append(sos)
        #Accepts mutation only if sum of squares value decreases
        if sos <= sos_list[-1]:
            output_file_name = "three_genes_replicated_" + str(i) + ".tsv"
            save_file(nf, output_file_name)
            save_file(genome_tracker, genome_output_file_name)
            sos_list.append(sos)
            gen+=1
    else:
        #Calculate fitness of new mutation
        probability = calc_fitness(f_new, f_old, Ne)
        #f_old
        genome_tracker['value'].values[12] = probability
        if probability > random.random():
            #Accepting mutation
            simulate_genome(genome_tracker, output_file_name)
            #Taking in new file and removing unnecessary rows and columns
            nf = pd.read_csv("three_genes_replicated.tsv", header=0, sep='\t')
            nf = edit_new_file(nf)
            sos = calc_sum_of_squares(df, nf)
            all_sos_list.append(sos)
            #Accepts mutation only if sum of squares value decreases
            if sos <= sos_list[-1]:
                output_file_name = "three_genes_replicated_" + str(i) + ".tsv"
                save_file(nf, output_file_name)
                save_file(genome_coordinates, genome_output_file_name)
                sos_list.append(sos)
                gen+=1

def alter_poly_strength(genome_tracker, which_promoter):

    poly_sigma = 1e5

    if which_promoter == "promoter":
        #Determining polymerase strength
        poly_eps = np.random.normal(mu, poly_sigma)
        #Determines new polymerase strength to reduce sum of squares value
        pol_strength = genome_tracker['new_strength'].values[3] + poly_eps
        while pol_strength < 0:
            poly_eps = np.random.normal(mu, poly_sigma)
            pol_strength = genome_tracker['new_strength'].values[3] + poly_eps
            if pol_strength > 3e10:
                poly_eps = np.random.normal(mu, poly_sigma)
                pol_strength = genome_tracker['new_strength'].values[3] + poly_eps
        genome_tracker['new_strength'].values[3] = pol_strength

    if which_promoter == "promoter1":
        #Determining polymerase strength
        poly_eps = np.random.normal(mu, poly_sigma)
        #Determines new polymerase strength to reduce sum of squares value
        pol1_strength = genome_tracker['new_strength'].values[4] + poly_eps
        while pol1_strength < 0:
            poly_eps = np.random.normal(mu, poly_sigma)
            pol1_strength = genome_tracker['new_strength'].values[4] + poly_eps
            if pol1_strength > 3e10:
                poly_eps = np.random.normal(mu, poly_sigma)
                pol1_strength = genome_tracker['new_strength'].values[4] + poly_eps
        genome_tracker['new_strength'].values[4] = pol1_strength

    if which_promoter == "promoter2":
        #Determining polymerase strength
        poly_eps = np.random.normal(mu, poly_sigma)
        #Determines new polymerase strength to reduce sum of squares value
        pol2_strength = genome_tracker['new_strength'].values[5] + poly_eps
        while pol2_strength < 0:
            poly_eps = np.random.normal(mu, poly_sigma)
            pol2_strength = genome_tracker['new_strength'].values[5] + poly_eps
            if pol2_strength > 3e10:
                poly_eps = np.random.normal(mu, poly_sigma)
                pol2_strength = genome_tracker['new_strength'].values[5] + poly_eps
        genome_tracker['new_strength'].values[5] = pol2_strength

def alter_term_efficiency(genome_tracker, which_terminator):

    term_sigma = 1.0

    if which_terminator == "terminator1":
        #Determining terminator polymerase efficiency rate
        term_eps = np.random.normal(mu, term_sigma)
        #Determines new terminator efficiency value to reduce sum of squares value
        term1_efficiency = genome_tracker['new_strength'].values[9] + term_eps
        while term1_efficiency < 0 or term1_efficiency > 1:
            term_eps = np.random.normal(mu, term_sigma)
            term1_efficiency = genome_tracker['new_strength'].values[9] + term_eps
        genome_tracker['new_strength'].values[9] = term1_efficiency

    if which_terminator == "terminator2":
        #Determining terminator polymerase efficiency rate
        term_eps = np.random.normal(mu, term_sigma)
        #Determines new terminator efficiency value to reduce sum of squares value
        term2_efficiency = genome_tracker['new_strength'].values[10] + term_eps
        while term2_efficiency < 0 or term2_efficiency > 1:
            term_eps = np.random.normal(mu, term_sigma)
            term2_efficiency = genome_tracker['new_strength'].values[10] + term_eps
        genome_tracker['new_strength'].values[10] = term2_efficiency

    if which_terminator == "terminator3":
        #Determining terminator polymerase efficiency rate
        term_eps = np.random.normal(mu, term_sigma)
        #Determines new terminator efficiency value to reduce sum of squares value
        term3_efficiency = genome_tracker['new_strength'].values[11] + term_eps
        while term3_efficiency < 0 or term3_efficiency > 1:
            term_eps = np.random.normal(mu, term_sigma)
            term3_efficiency = genome_tracker['new_strength'].values[11] + term_eps
        genome_tracker['new_strength'].values[11] = term3_efficiency

def add_promoter(genome_tracker):

    promoter_possibilities = ["promoter1", "promoter2"]
    chosen_promoter = random.choice(promoter_possibilities)
    promoter1_slots = ['A', 'B']
    promoter2_slots = ['A', 'B']

    if chosen_promoter == "promoter1":
        #Adding in a promoter between genes 1 and 2
        if (genome_tracker['start'].values[6] >= genome_tracker['start'].values[15]) and (genome_tracker['start'].values[6] <= genome_tracker['stop'].values[15]):
            promoter1_slots.remove('A')
        if (genome_tracker['start'].values[9] >= genome_tracker['start'].values[15]) and (genome_tracker['start'].values[9] <= genome_tracker['stop'].values[15]):
            promoter1_slots.remove('A')
        if (genome_tracker['start'].values[6] >= genome_tracker['start'].values[16]) and (genome_tracker['start'].values[6] <= genome_tracker['stop'].values[16]):
            promoter1_slots.remove('B')
        if (genome_tracker['start'].values[9] >= genome_tracker['start'].values[16]) and (genome_tracker['stop'].values[9] <= genome_tracker['stop'].values[16]):
            promoter1_slots.remove('B')
        if promoter1_slots == []:
            genome_tracker['start'].values[4] = 0
            genome_tracker['stop'].values[4] = 0
        else:
            available_slot = random.choice(promoter1_slots)
            if available_slot == 'A':
                prom1_start = genome_tracker['start'].values[4] = random.randint(genome_tracker['start'].values[15], 123)
                prom1_stop = genome_tracker['stop'].values[4] = prom1_start + 9
            if available_slot == 'B':
                prom1_start = genome_tracker['start'].values[4] = random.randint(genome_tracker['start'].values[16], 134)
                prom1_stop = genome_tracker['stop'].values[4] = prom1_start + 9

        genome_tracker['previous_strength'].values[4] = genome_tracker['new_strength'].values[4]
        genome_tracker['new_strength'].values[4] = alter_poly_strength(genome_tracker, chosen_promoter)

    if chosen_promoter == "promoter2":
        #Adding promoter between genes 2 and 3
        if (genome_tracker['start'].values[7] >= genome_tracker['start'].values[18]) and (genome_tracker['start'].values[7] <= genome_tracker['stop'].values[18]):
            promoter2_slots.remove('A')
        if (genome_tracker['start'].values[10] >= genome_tracker['start'].values[18]) and (genome_tracker['start'].values[10] <= genome_tracker['stop'].values[18]):
            promoter2_slots.remove('A')
        if (genome_tracker['start'].values[7] >= genome_tracker['start'].values[19]) and (genome_tracker['start'].values[7] <= genome_tracker['stop'].values[19]):
            promoter2_slots.remove('B')
        if (genome_tracker['start'].values[10] >= genome_tracker['start'].values[19]) and (genome_tracker['start'].values[10] <= genome_tracker['stop'].values[19]):
            promoter2_slots.remove('B')
        if promoter2_slots == []:
            genome_tracker['start'].values[5] = 0
            genome_tracker['stop'].values[5] = 0
        else:
            available_slot = random.choice(promoter2_slots)
            if available_slot == 'A':
                prom2_start = genome_tracker['start'].values[5] = random.randint(genome_tracker['start'].values[18], 282)
                prom2_stop = genome_tracker['stop'].values[5] = prom2_start + 9
            if available_slot == 'B':
                prom2_start = genome_tracker['start'].values[5] = random.randint(genome_tracker['start'].values[19], 293)
                prom2_stop = genome_tracker['stop'].values[5] = prom2_start + 9

        genome_tracker['previous_strength'].values[5] = genome_tracker['new_strength'].values[5]
        genome_tracker['new_strength'].values[5] = alter_poly_strength(genome_tracker, chosen_promoter)

def remove_promoter(genome_tracker):

    promoter_possibilities = ["promoter1", "promoter2"]
    chosen_promoter = random.choice(promoter_possibilities)

    if chosen_promoter == "promoter1":
        #Removing promoter between genes 1 and 2
        genome_tracker['start'].values[4] = 0
        genome_tracker['stop'].values[4] = 0

    if chosen_promoter == "promoter2":
        #Removing promoter between genes 2 and 3
        genome_tracker['start'].values[5] = 0
        genome_tracker['stop'].values[5] = 0

def add_rnase(genome_tracker):

    rnase_possibilities = ["rnase1", "rnase2", "rnase3"]
    chosen_rnase = random.choice(rnase_possibilities)
    rnase_slots = ['A', 'B', 'C']

    if chosen_rnase == "rnase1":
        #Adds rnase after first promoter
        rnase1_start = genome_tracker['start'].values[6] = random.randint(genome_tracker['start'].values[14], 15)
        rnase1_stop = genome_tracker['stop'].values[6] = rnase1_start + 10

    if chosen_rnase == "rnase2":
        #Adds rnase after second promoter
        if (genome_tracker['start'].values[4] >= genome_tracker['start'].values[15]) and (genome_tracker['start'].values[4] <= genome_tracker['stop'].values[15]):
            rnase_slots.remove('A')
        if (genome_tracker['start'].values[9] >= genome_tracker['start'].values[15]) and (genome_tracker['start'].values[9] <= genome_tracker['stop'].values[15]):
            rnase_slots.remove('A')
        if (genome_tracker['start'].values[4] >= genome_tracker['start'].values[16]) and (genome_tracker['start'].values[4] <= genome_tracker['stop'].values[16]):
            rnase_slots.remove('B')
        if (genome_tracker['start'].values[9] >= genome_tracker['start'].values[16]) and (genome_tracker['start'].values[9] <= genome_tracker['stop'].values[16]):
            rnase_slots.remove('B')
        if (genome_tracker['start'].values[9] >= genome_tracker['start'].values[17]) and (genome_tracker['start'].values[9] <= genome_tracker['stop'].values[17]):
            rnase_slots.remove('C')
        if rnase_slots == []:
            genome_tracker['start'].values[7] = 0
            genome_tracker['stop'].values[7] = 0
        else:
            available_slot = random.choice(rnase_slots)
            if available_slot == 'A':
                rnase2_start = genome_tracker['start'].values[7] = random.randint(genome_tracker['start'].values[15], 122)
                rnase2_stop = genome_tracker['stop'].values[7] = rnase2_start + 10
            if available_slot == 'B':
                rnase2_start = genome_tracker['start'].values[7] = random.randint(genome_tracker['start'].values[16], 133)
                rnase2_stop = genome_tracker['stop'].values[7] = rnase2_start + 10
            if available_slot == 'C':
                rnase2_start = genome_tracker['start'].values[7] = random.randint(genome_tracker['start'].values[17], 149)
                rnase2_stop = genome_tracker['stop'].values[7] = rnase2_start + 10

        if chosen_rnase == "rnase3":
            #Adds rnase after third promoter
            if (genome_tracker['start'].values[5] >= genome_tracker['start'].values[18]) and (genome_tracker['start'].values[5] <= genome_tracker['stop'].values[18]):
                rnase_slots.remove('A')
            if (genome_tracker['start'].values[10] >= genome_tracker['start'].values[18]) and (genome_tracker['start'].values[10] <= genome_tracker['stop'].values[18]):
                rnase_slots.remove('A')
            if (genome_tracker['start'].values[5] >= genome_tracker['start'].values[19]) and (genome_tracker['start'].values[5] <= genome_tracker['stop'].values[19]):
                rnase_slots.remove('B')
            if (genome_tracker['start'].values[10] >= genome_tracker['start'].values[19]) and (genome_tracker['start'].values[10] <= genome_tracker['stop'].values[19]):
                rnase_slots.remove('B')
            if (genome_tracker['start'].values[10] >= genome_tracker['start'].values[20]) and (genome_tracker['start'].values[10] <= genome_tracker['stop'].values[20]):
                rnase_slots.remove('C')
            if rnase_slots == []:
                genome_tracker['start'].values[8] = 0
                genome_tracker['stop'].values[8] = 0
            else:
                available_slot = random.choice(rnase_slots)
                if available_slot == 'A':
                    rnase3_start = genome_tracker['start'].values[8] = random.randint(region2a_start, 281)
                    rnase3_stop = genome_tracker['stop'].values[8] = rnase3_start + 10
                if available_slot == 'B':
                    rnase3_start = genome_tracker['start'].values[8] = random.randint(region2b_start, 292)
                    rnase3_stop = genome_tracker['stop'].values[8] = rnase3_start + 10
                if available_slot == 'C':
                    rnase3_start = genome_tracker['start'].values[8] = random.randint(region2c_start, 308)
                    rnase3_stop = genome_tracker['stop'].values[8] = rnase3_start + 10

def remove_rnase(genome_tracker):

    rnase_possibilities = ["rnase1", "rnase2", "rnase3"]
    chosen_rnase = random.choice(rnase_possibilities)

    if chosen_rnase == "rnase1":
        #Removes rnase after first promoter
        genome_tracker['start'].values[6] = 0
        genome_tracker['stop'].values[6] = 0

    if chosen_rnase == "rnase2":
        #Removes rnase after second promoter
        genome_tracker['start'].values[7] = 0
        genome_tracker['stop'].values[7] = 0

    if chosen_rnase == "rnase3":
        #Removes rnase after third promoter
        genome_tracker['start'].values[8] = 0
        genome_tracker['stop'].values[8] = 0

def add_terminator(genome_tracker):

    terminator_possibilities = ["terminator1", "terminator2", "terminator3"]
    chosen_terminator = random.choice(terminator_possibilities)
    terminator_slots = ['A', 'B', 'C']

    if chosen_terminator == "terminator1":
        #Adds terminator after first gene
        if (genome_tracker['start'].values[4] >= genome_tracker['start'].values[15]) and (genome_tracker['start'].values[4] <= genome_tracker['stop'].values[15]):
            terminator_slots.remove('A')
        if (genome_tracker['start'].values[6] >= genome_tracker['start'].values[15]) and (genome_tracker['start'].values[6] <= genome_tracker['stop'].values[15]):
            terminator_slots.remove('A')
        if (genome_tracker['start'].values[4] >= genome_tracker['start'].values[16]) and (genome_tracker['start'].values[4] <= genome_tracker['stop'].values[16]):
            terminator_slots.remove('B')
        if (genome_tracker['start'].values[6] >= genome_tracker['start'].values[16]) and (genome_tracker['start'].values[6] <= genome_tracker['stop'].values[16]):
            terminator_slots.remove('B')
        if (genome_tracker['start'].values[6] >= genome_tracker['start'].values[17]) and (genome_tracker['start'].values[6] <= genome_tracker['stop'].values[17]):
            terminator_slots.remove('C')
        if terminator_slots == []:
            genome_tracker['start'].values[9] = 0
            genome_tracker['stop'].values[9] = 0
        else:
            available_slot = random.choice(terminator_slots)
            if available_slot == 'A':
                term1_start = genome_tracker['start'].values[9] = random.randint(genome_tracker['start'].values[15], 131)
                term1_stop = genome_tracker['stop'].values[9] = term1_start + 1
            if available_slot == 'B':
                term1_start = genome_tracker['start'].values[9] = random.randint(genome_tracker['start'].values[16], 142)
                term1_stop = genome_tracker['stop'].values[9] = term1_start + 1
            if available_slot == 'C':
                term1_start = genome_tracker['start'].values[9] = random.randint(genome_tracker['start'].values[17], 158)
                term1_stop = genome_tracker['stop'].values[9] = term1_start + 1

        genome_tracker['previous_strength'].values[9] = genome_tracker['new_strength'].values[9]
        genome_tracker['new_strength'].values[9] = alter_term_efficiency(genome_tracker, chosen_terminator)

    if chosen_terminator == "terminator2":
        #Adds terminator after second gene
        if (genome_tracker['start'].values[5] >= genome_tracker['start'].values[18]) and (genome_tracker['start'].values[5] <= genome_tracker['stop'].values[18]):
            terminator_slots.remove('A')
        if (genome_tracker['start'].values[7] >= genome_tracker['start'].values[18]) and (genome_tracker['start'].values[7] <= genome_tracker['stop'].values[18]):
            terminator_slots.remove('A')
        if (genome_tracker['start'].values[5] >= genome_tracker['start'].values[19]) and (genome_tracker['start'].values[5] <= genome_tracker['stop'].values[19]):
            terminator_slots.remove('B')
        if (genome_tracker['start'].values[7] >= genome_tracker['start'].values[19]) and (genome_tracker['start'].values[7] <= genome_tracker['stop'].values[19]):
            terminator_slots.remove('B')
        if (genome_tracker['start'].values[7] >= genome_tracker['start'].values[20]) and (genome_tracker['start'].values[7] <= genome_tracker['stop'].values[20]):
            terminator_slots.remove('C')
        if terminator_slots == []:
            genome_tracker['start'].values[10] = 0
            genome_tracker['stop'].values[10] = 0
        else:
            available_slot = random.choice(terminator_slots)
            if available_slot == 'A':
                term2_start = genome_tracker['start'].values[10] = random.randint(genome_tracker['start'].values[15], 290)
                term2_stop = genome_tracker['stop'].values[10] = term2_start + 1
            if available_slot == 'B':
                term2_start = genome_tracker['start'].values[10] = random.randint(genome_tracker['start'].values[16], 301)
                term2_stop = genome_tracker['stop'].values[10] = term2_start + 1
            if available_slot == 'C':
                term2_start = genome_tracker['start'].values[10] = random.randint(genome_tracker['start'].values[17], 317)
                term2_stop = genome_tracker['stop'].values[10] = term2_start + 1

        genome_tracker['previous_strength'].values[10] = genome_tracker['new_strength'].values[10]
        genome_tracker['new_strength'].values[10] = alter_term_efficiency(genome_tracker, chosen_terminator)

    if chosen_terminator == "terminator3":
        #Adds terminator after third gene
        genome_tracker['start'].values[11] = 449
        genome_tracker['stop'].values[11] = 450

        genome_tracker['previous_strength'].values[11] = genome_tracker['new_strength'].values[11]
        genome_tracker['new_strength'].values[11] = alter_term_efficiency(genome_tracker, chosen_terminator)

def remove_terminator(genome_tracker):

    terminator_possibilities = ["terminator1", "terminator2", "terminator3"]
    chosen_terminator = random.choice(terminator_possibilities)

    if chosen_terminator == "terminator1":
        #Removes terminator after first gene
        genome_tracker['start'].values[9] = 0
        genome_tracker['stop'].values[9] = 0
        genome_tracker['previous_strength'].values[9] = genome_tracker['new_strength'].values[9]
        genome_tracker['new_strength'].values[9] = 0.0

    if chosen_terminator == "terminator2":
        #Removes terminator after second gene
        genome_tracker['start'].values[10] = 0
        genome_tracker['stop'].values[10] = 0
        genome_tracker['previous_strength'].values[10] = genome_tracker['new_strength'].values[10]
        genome_tracker['new_strength'].values[10] = 0.0

    if chosen_terminator == "terminator3":
        #Removes terminator after third gene
        genome_tracker['start'].values[11] = 0
        genome_tracker['stop'].values[11] = 0
        genome_tracker['previous_strength'].values[11] = genome_tracker['new_strength'].values[11]
        genome_tracker['new_strength'].values[11] = 0.0

def save_file(file, file_name):

    f = pd.DataFrame(file)
    f.to_csv(file_name, sep='\t')

def simulate_genome(genome_tracker, file_name):

    genome_tracker['start'] = genome_tracker['start'].astype(int)
    genome_tracker['stop'] = genome_tracker['stop'].astype(int)
    sim = pt.Model(cell_volume=8e-16)
    sim.seed(34)
    sim.add_polymerase(name="rnapol", copy_number=4, speed=40, footprint=10)
    sim.add_ribosome(copy_number=100, speed=30, footprint=10)

    plasmid = pt.Genome(name="plasmid", length=450,
                        transcript_degradation_rate=1e-2,
                        transcript_degradation_rate_ext=1e-2,
                        rnase_speed=20,
                        rnase_footprint=10)
    plasmid.add_promoter(name="p1", start=1, stop=10,
                         interactions={"rnapol": genome_tracker['new_strength'].values[3]})
    if genome_tracker['start'].values[4] > 0:
        plasmid.add_promoter(name="p2", start=genome_tracker['start'].values[4], stop=genome_tracker['stop'].values[4],
                             interactions={"rnapol": genome_tracker['new_strength'].values[4]})
    if genome_tracker['start'].values[5] > 0:
        plasmid.add_promoter(name="p3", start=genome_tracker['start'].values[5], stop=genome_tracker['stop'].values[5],
                             interactions={"rnapol": genome_tracker['new_strength'].values[5]})
    if genome_tracker['start'].values[6] > 0:
        plasmid.add_rnase_site(start=genome_tracker['start'].values[6], stop=genome_tracker['stop'].values[6])
    if genome_tracker['start'].values[7] > 0:
        plasmid.add_rnase_site(start=genome_tracker['start'].values[7], stop=genome_tracker['stop'].values[7])
    if genome_tracker['start'].values[8] > 0:
        plasmid.add_rnase_site(start=genome_tracker['start'].values[8], stop=genome_tracker['stop'].values[8])
    if genome_tracker['start'].values[9] > 0:
        plasmid.add_terminator(name="t1", start=genome_tracker['start'].values[9], stop=genome_tracker['stop'].values[9],
                               efficiency={"rnapol": genome_tracker['new_strength'].values[9]})
    if genome_tracker['start'].values[10] > 0:
        plasmid.add_terminator(name="t2", start=genome_tracker['start'].values[10], stop=genome_tracker['stop'].values[10],
                               efficiency={"rnapol": genome_tracker['new_strength'].values[10]})
    if genome_tracker['start'].values[11] > 0:
        plasmid.add_terminator(name="t3", start=genome_tracker['start'].values[11], stop=genome_tracker['stop'].values[11],
                               efficiency={"rnapol": genome_tracker['new_strength'].values[11]})
    plasmid.add_gene(name="proteinX", start=genome_tracker['start'].values[0], stop=genome_tracker['stop'].values[0],
                     rbs_start=(genome_tracker['start'].values[0]-15), rbs_stop=genome_tracker['start'].values[0], rbs_strength=1e7)

    plasmid.add_gene(name="proteinY", start=genome_tracker['start'].values[1], stop=genome_tracker['stop'].values[1],
                     rbs_start=(genome_tracker['start'].values[1]-15), rbs_stop=genome_tracker['start'].values[1], rbs_strength=1e7)
    plasmid.add_gene(name="proteinZ", start=genome_tracker['start'].values[2], stop=genome_tracker['stop'].values[2],
                     rbs_start=(genome_tracker['start'].values[2]-15), rbs_stop=genome_tracker['start'].values[2], rbs_strength=1e7)
    sim.register_genome(plasmid)
    sim.simulate(time_limit=240, time_step=1, output=file_name)

if __name__ == '__main__':
    main()
