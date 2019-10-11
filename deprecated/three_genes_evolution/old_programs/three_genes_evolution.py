import pinetree as pt
import pandas as pd
import random
import numpy as np
import string
import filecmp
import math
import os
import itertools

#Need to add in the following functions: expand genome when adding in components, simulated annealing

def main():

    global Ne
    global gen
    global mu
    global sigma
    sos_list = [60000]
    poly_list = []
    poly1_list = []
    poly2_list = []
    all_poly_list = []
    all_sos_list = []
    term1_efficiency_list = []
    term2_efficiency_list = []
    term3_efficiency_list = []
    all_term_list = []
    ribo1_strength_list = []
    ribo2_strength_list = []
    ribo3_strength_list = []
    gene1_start_list = []
    gene2_start_list = []
    gene3_start_list = []
    gene1_stop_list = []
    gene2_stop_list = []
    gene3_stop_list = []
    gene1_list = [0, 0]
    gene2_list = [0, 0]
    gene3_list = [0, 0]
    term1_start_list = [0]
    term1_stop_list = [0]
    term2_start_list = [0]
    term2_stop_list = [0]
    term3_start_list = [0]
    term3_stop_list = [0]
    term1_list = [0, 0]
    term2_list = [0, 0]
    term3_list = [0, 0]
    rnase1_start_list = [0]
    rnase1_stop_list = [0]
    rnase2_start_list = [0]
    rnase2_stop_list = [0]
    rnase3_start_list = [0]
    rnase3_stop_list = [0]
    rnase1_list = [0, 0]
    rnase2_list = [0, 0]
    rnase3_list = [0, 0]
    prom1_start_list = [0]
    prom1_stop_list = [0]
    prom2_start_list = [0]
    prom2_stop_list = [0]
    prom1_list = [0, 0]
    prom2_list = [0, 0]
    iter_list = []
    f_old = 1.0
    mu = 0.0
    sigma = 1.0
    Ne = 10
    i = 0
    gen = 0
    new_gene1_start = 26
    new_gene1_stop = 121
    new_gene2_start = 159
    new_gene2_stop = 280
    new_gene3_start = 319
    new_gene3_stop = 449
    new_prom1_start = 0
    new_prom1_stop = 0
    new_prom2_start = 0
    new_prom2_stop = 0
    new_rnase1_start = 0
    new_rnase1_stop = 0
    new_rnase2_start = 0
    new_rnase2_stop = 0
    new_rnase3_start = 0
    new_rnase3_stop = 0
    new_term1_start = 0
    new_term1_stop = 0
    new_term2_start = 0
    new_term2_stop = 0
    new_term3_start = 0
    new_term3_stop = 0

    #Taking in target tsv file
    df_name = input("Enter tsv file name: ")
    df = pd.read_csv(df_name, header=0, sep='\t')
    df = edit_target_file(df, df_name)

    #Evolution program
    new_pol_strength = random.randint(0, 3e10)
    new_pol1_strength = random.randint(0, 3e10)
    new_pol2_strength = random.randint(0, 3e10)
    new_term1_efficiency = random.uniform(0.0, 1.0)
    new_term2_efficiency = random.uniform(0.0, 1.0)
    new_term3_efficiency = random.uniform(0.0, 1.0)
    poly_list.append(new_pol_strength)
    poly1_list.append(new_pol1_strength)
    poly2_list.append(new_pol2_strength)
    term1_efficiency_list.append(new_term1_efficiency)
    term2_efficiency_list.append(new_term2_efficiency)
    term3_efficiency_list.append(new_term3_efficiency)

    possibilities = ["alter polymerase 1 strength", "add promoter", "remove promoter", "add rnase",
                      "remove rnase", "add terminator", "remove terminator"]

    while i < 5000:

        mutation = random.choice(possibilities)

        if mutation == "alter polymerase 1 strength":
            promoter = "promoter"
            new_pol_strength = alter_poly_strength(poly_list, poly1_list, poly2_list, promoter)

        if mutation == "alter gene length":
            gene_length = alter_gene_length(new_gene1_start, new_gene1_stop, new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop)
            new_gene1_start = gene1_start_list[-1]
            new_gene1_stop = gene1_stop_list[-1]
            new_gene2_start = gene2_start_list[-1]
            new_gene2_stop = gene2_stop_list[-1]
            new_gene3_start = gene3_start_list[-1]
            new_gene3_stop = gene3_stop_list[-1]

        if mutation == "add promoter":
            new_promoter = add_promoter(new_prom1_start, new_prom1_stop, new_prom2_start, new_prom2_stop, prom1_start_list, prom1_stop_list,
                                         prom2_start_list, prom2_stop_list, term1_start_list, term1_stop_list, term2_start_list, term2_stop_list,
                                         rnase1_start_list, rnase1_stop_list, rnase2_start_list, rnase2_stop_list, poly_list, poly1_list, poly2_list)
            new_prom1_start = prom1_start_list[-1]
            new_prom1_stop = prom1_stop_list[-1]
            new_prom2_start = prom2_start_list[-1]
            new_prom2_stop = prom2_stop_list[-1]
            new_pol1_strength = new_promoter[0]
            new_pol2_strength = new_promoter[1]

        if mutation == "remove promoter":
            old_promoter = remove_promoter(new_prom1_start, new_prom1_stop, new_prom2_start, new_prom2_stop, prom1_start_list, prom1_stop_list,
                                            prom2_start_list, prom2_stop_list)
            new_prom1_start = prom1_start_list[-1]
            new_prom1_stop = prom1_stop_list[-1]
            new_prom2_start = prom2_start_list[-1]
            new_prom2_stop = prom2_stop_list[-1]

        if mutation == "add rnase":
            new_rnase = add_rnase(new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop, new_rnase3_start, new_rnase3_stop,
                                    prom1_start_list, prom1_stop_list, prom2_start_list, prom2_stop_list, term1_start_list, term1_stop_list,
                                    term2_start_list, term2_stop_list, rnase1_start_list, rnase1_stop_list, rnase2_start_list, rnase2_stop_list,
                                    rnase3_start_list, rnase3_stop_list)
            new_rnase1_start = rnase1_start_list[-1]
            new_rnase1_stop = rnase1_stop_list[-1]
            new_rnase2_start = rnase2_start_list[-1]
            new_rnase2_stop = rnase2_stop_list[-1]
            new_rnase3_start = rnase3_start_list[-1]
            new_rnase3_stop = rnase3_stop_list[-1]

        if mutation == "remove rnase":
            old_rnase = remove_rnase(new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop, new_rnase3_start, new_rnase3_stop,
                                        rnase1_start_list, rnase1_stop_list, rnase2_start_list, rnase2_stop_list, rnase3_start_list, rnase3_stop_list)
            new_rnase1_start = rnase1_start_list[-1]
            new_rnase1_stop = rnase1_stop_list[-1]
            new_rnase2_start = rnase2_start_list[-1]
            new_rnase2_stop = rnase2_stop_list[-1]
            new_rnase3_start = rnase3_start_list[-1]
            new_rnase3_stop = rnase3_stop_list[-1]

        if mutation == "add terminator":
            new_terminator = add_terminator(new_term1_start, new_term1_stop, new_term2_start, new_term2_stop, new_term3_start, new_term3_stop,
                                                prom1_start_list, prom1_stop_list, prom2_start_list, prom2_stop_list, term1_start_list, term1_stop_list,
                                                term2_start_list, term2_stop_list, term3_start_list, term3_stop_list, rnase1_start_list, rnase1_stop_list,
                                                rnase2_start_list, rnase2_stop_list, term1_efficiency_list, term2_efficiency_list, term3_efficiency_list)
            new_term1_start = term1_start_list[-1]
            new_term1_stop = term1_stop_list[-1]
            new_term2_start = term2_start_list[-1]
            new_term2_stop = term2_stop_list[-1]
            new_term3_start = term3_start_list[-1]
            new_term3_stop = term3_stop_list[-1]
            new_term1_efficiency = new_terminator[0]
            new_term2_efficiency = new_terminator[1]
            new_term3_efficiency = new_terminator[2]

        if mutation == "remove terminator":
            old_terminator = remove_terminator(new_term1_start, new_term1_stop, new_term2_start, new_term2_stop, new_term3_start, new_term3_stop,
                                                    term1_start_list, term1_stop_list, term2_start_list, term2_stop_list, term3_start_list, term3_stop_list)
            new_term1_start = term1_start_list[-1]
            new_term1_stop = term1_stop_list[-1]
            new_term2_start = term2_start_list[-1]
            new_term2_stop = term2_stop_list[-1]
            new_term3_start = term3_start_list[-1]
            new_term3_stop = term3_stop_list[-1]

        eps = np.random.normal(mu, sigma)
        f_new = f_old * (1.0 + eps)
        test = accept_mutation(df, f_old, f_new, new_pol_strength, new_pol1_strength, new_pol2_strength, new_term1_efficiency, new_term2_efficiency,
                                    new_term3_efficiency, new_gene1_start, new_gene1_stop, new_gene2_start,
                                    new_gene2_stop, new_gene3_start, new_gene3_stop, new_prom1_start, new_prom1_stop, new_prom2_start,
                                    new_prom2_stop, new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop, new_rnase3_start,
                                    new_rnase3_stop, new_term1_start, new_term1_stop, new_term2_start, new_term2_stop, new_term3_start,
                                    new_term3_stop, all_sos_list, sos_list, poly_list, poly1_list, poly2_list, term1_efficiency_list, term2_efficiency_list,
                                    term3_efficiency_list, gene1_list, gene2_list, gene3_list, prom1_list, prom2_list, term1_list, term2_list, term3_list,
                                    rnase1_list, rnase2_list, rnase3_list, i, iter_list)
        if test > 0:
            f_old = test
        '''if sos_list[-1] == 0:
            break'''

        if i == 0:
            output_file = "gen_0_data.tsv"
            three_genome.recreated_genome(new_pol_strength, new_pol1_strength, new_pol2_strength, new_term1_efficiency,
                                            new_term2_efficiency, new_term3_efficiency, new_gene1_start, new_gene1_stop,
                                            new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop,
                                            new_prom1_start, new_prom1_stop, new_prom2_start, new_prom2_stop,
                                            new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop,
                                            new_rnase3_start, new_rnase3_stop, new_term1_start, new_term1_stop,
                                            new_term2_start, new_term2_stop, new_term3_start, new_term3_stop, output_file)
            #Gets genome coordinates at the 1st iteration
            iteration = "0"
            gene1_list[0] = new_gene1_start
            gene1_list[1] = new_gene1_stop
            gene2_list[0] = new_gene2_start
            gene2_list[1] = new_gene2_stop
            gene3_list[0] = new_gene3_start
            gene3_list[1] = new_gene3_stop
            prom1_list[0] = new_prom1_start
            prom1_list[1] = new_prom1_stop
            prom2_list[0] = new_prom2_start
            prom2_list[1] = new_prom2_stop
            term1_list[0] = new_term1_start
            term1_list[1] = new_term1_stop
            term2_list[0] = new_term2_start
            term2_list[1] = new_term2_stop
            term3_list[0] = new_term3_start
            term3_list[1] = new_term3_stop
            rnase1_list[0] = new_rnase1_start
            rnase1_list[1] = new_rnase1_stop
            rnase2_list[0] = new_rnase2_start
            rnase2_list[1] = new_rnase2_stop
            rnase3_list[0] = new_rnase3_start
            rnase3_list[1] = new_rnase3_stop
            get_coordinates(gene1_list, gene2_list, gene3_list, prom1_list, prom2_list, term1_list, term2_list, term3_list,
                                rnase1_list, rnase2_list, rnase3_list, poly_list, poly1_list, poly2_list, term1_efficiency_list,
                                term2_efficiency_list, term3_efficiency_list, iteration)

        if i == 4999:
            output_file = "gen_5000_data.tsv"
            three_genome.recreated_genome(new_pol_strength, new_pol1_strength, new_pol2_strength, new_term1_efficiency,
                                            new_term2_efficiency, new_term3_efficiency, new_gene1_start, new_gene1_stop,
                                            new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop,
                                            new_prom1_start, new_prom1_stop, new_prom2_start, new_prom2_stop,
                                            new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop,
                                            new_rnase3_start, new_rnase3_stop, new_term1_start, new_term1_stop,
                                            new_term2_start, new_term2_stop, new_term3_start, new_term3_stop, output_file)
            #Gets genome coordinates at the 500th iteration
            iteration = "5000"
            gene1_list[0] = new_gene1_start
            gene1_list[1] = new_gene1_stop
            gene2_list[0] = new_gene2_start
            gene2_list[1] = new_gene2_stop
            gene3_list[0] = new_gene3_start
            gene3_list[1] = new_gene3_stop
            prom1_list[0] = new_prom1_start
            prom1_list[1] = new_prom1_stop
            prom2_list[0] = new_prom2_start
            prom2_list[1] = new_prom2_stop
            term1_list[0] = new_term1_start
            term1_list[1] = new_term1_stop
            term2_list[0] = new_term2_start
            term2_list[1] = new_term2_stop
            term3_list[0] = new_term3_start
            term3_list[1] = new_term3_stop
            rnase1_list[0] = new_rnase1_start
            rnase1_list[1] = new_rnase1_stop
            rnase2_list[0] = new_rnase2_start
            rnase2_list[1] = new_rnase2_stop
            rnase3_list[0] = new_rnase3_start
            rnase3_list[1] = new_rnase3_stop
            get_coordinates(gene1_list, gene2_list, gene3_list, prom1_list, prom2_list, term1_list, term2_list, term3_list,
                                rnase1_list, rnase2_list, rnase3_list, poly_list, poly1_list, poly2_list, term1_efficiency_list,
                                term2_efficiency_list, term3_efficiency_list, iteration)

        if i == 9999:
            output_file = "gen_10000_data.tsv"
            three_genome.recreated_genome(new_pol_strength, new_pol1_strength, new_pol2_strength, new_term1_efficiency,
                                            new_term2_efficiency, new_term3_efficiency, new_gene1_start, new_gene1_stop,
                                            new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop,
                                            new_prom1_start, new_prom1_stop, new_prom2_start, new_prom2_stop,
                                            new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop,
                                            new_rnase3_start, new_rnase3_stop, new_term1_start, new_term1_stop,
                                            new_term2_start, new_term2_stop, new_term3_start, new_term3_stop, output_file)
            #Gets genome coordinates at the 1000th iteration
            iteration = "10000"
            gene1_list[0] = new_gene1_start
            gene1_list[1] = new_gene1_stop
            gene2_list[0] = new_gene2_start
            gene2_list[1] = new_gene2_stop
            gene3_list[0] = new_gene3_start
            gene3_list[1] = new_gene3_stop
            prom1_list[0] = new_prom1_start
            prom1_list[1] = new_prom1_stop
            prom2_list[0] = new_prom2_start
            prom2_list[1] = new_prom2_stop
            term1_list[0] = new_term1_start
            term1_list[1] = new_term1_stop
            term2_list[0] = new_term2_start
            term2_list[1] = new_term2_stop
            term3_list[0] = new_term3_start
            term3_list[1] = new_term3_stop
            rnase1_list[0] = new_rnase1_start
            rnase1_list[1] = new_rnase1_stop
            rnase2_list[0] = new_rnase2_start
            rnase2_list[1] = new_rnase2_stop
            rnase3_list[0] = new_rnase3_start
            rnase3_list[1] = new_rnase3_stop
            get_coordinates(gene1_list, gene2_list, gene3_list, prom1_list, prom2_list, term1_list, term2_list, term3_list,
                                rnase1_list, rnase2_list, rnase3_list, poly_list, poly1_list, poly2_list, term1_efficiency_list,
                                term2_efficiency_list, term3_efficiency_list, iteration)

        i+=1
        print("i = ", i)


    #Exported tsv files
    '''iteration+=1
    sim_number = iteration
    path = "~/pinetree-toys/three_genes_evolution/" + string(sim_number)
    if not os.path.exists(path):
        os.makedirs(path)'''
    all_sos_dataframe = pd.DataFrame(all_sos_list, columns=["Sum_of_Squares"])
    export_csv = all_sos_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/all_sos_data.tsv", index=False)
    sos_list.remove(60000)
    sos_dataframe = pd.DataFrame(data=sos_list, columns=["Sum_of_Squares"])
    export_csv = sos_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/sos_data.tsv", index=False)
    iter_dataframe = pd.DataFrame(data=iter_list, columns=["Iteration"])
    export_csv = iter_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/sos_iter_data.tsv", index=False)
    all_poly_dataframe = pd.DataFrame(all_poly_list, columns=["Polymerase_Rate"])
    export_csv = all_poly_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/all_poly_data.tsv", index=False)
    poly_dataframe = pd.DataFrame(poly_list, columns=["Polymerase_Rate"])
    export_csv = poly_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/poly_data.tsv", index=False)
    poly1_dataframe = pd.DataFrame(poly1_list, columns=["Polymerase_Rate"])
    export_csv = poly1_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/poly1_data.tsv", index=False)
    poly2_dataframe = pd.DataFrame(poly2_list, columns=["Polymerase_Rate"])
    export_csv = poly2_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/poly2_data.tsv", index=False)
    all_term_dataframe = pd.DataFrame(all_term_list, columns=["Terminator Efficiency Rate"])
    export_csv = all_term_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/all_term_data.tsv", index=False)
    term1_efficiency_dataframe = pd.DataFrame(term1_efficiency_list, columns=["Terminator1_Efficiency_Rate"])
    export_csv = term1_efficiency_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/term1_efficiency_data.tsv", index=False)
    term2_efficiency_dataframe = pd.DataFrame(term2_efficiency_list, columns=["Terminator2_Efficiency_Rate"])
    export_csv = term2_efficiency_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/term2_efficiency_data.tsv", index=False)
    term3_efficiency_dataframe = pd.DataFrame(term3_efficiency_list, columns=["Terminator3_Efficiency_Rate"])
    export_csv = term3_efficiency_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/term3_efficiency_data.tsv", index=False)
    gene1_start_dataframe = pd.DataFrame(gene1_start_list, columns=["Gene_Start"])
    export_csv = gene1_start_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/gene1_start_data.tsv", index=False)
    gene2_start_dataframe = pd.DataFrame(gene2_start_list, columns=["Gene_Start"])
    export_csv = gene2_start_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/gene2_start_data.tsv", index=False)
    gene3_start_dataframe = pd.DataFrame(gene3_start_list, columns=["Gene Start"])
    export_csv = gene3_start_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/gene3_start_data.tsv", index=False)
    gene1_stop_dataframe = pd.DataFrame(gene1_stop_list, columns=["Gene_Stop"])
    export_csv = gene1_stop_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/gene1_stop_data.tsv", index=False)
    gene2_stop_dataframe = pd.DataFrame(gene2_stop_list, columns=["Gene_Stop"])
    export_csv = gene2_stop_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/gene2_stop_data.tsv", index=False)
    gene3_stop_dataframe = pd.DataFrame(gene3_stop_list, columns=["Gene_Stop"])
    export_csv = gene3_stop_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/gene3_stop_data.tsv", index=False)
    print("Generations = ", gen)

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
def sum_of_squares(target_file, new_file):

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

def accept_mutation(df, f_old, f_new, new_pol_strength, new_pol1_strength, new_pol2_strength, new_term1_efficiency, new_term2_efficiency,
                    new_term3_efficiency, new_gene1_start, new_gene1_stop, new_gene2_start,
                    new_gene2_stop, new_gene3_start, new_gene3_stop, promoter1_start, promoter1_stop, promoter2_start,
                    promoter2_stop, new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop, new_rnase3_start,
                    new_rnase3_stop, new_term1_start, new_term1_stop, new_term2_start, new_term2_stop, new_term3_start,
                    new_term3_stop, all_sos_list, sos_list, poly_list, poly1_list, poly2_list, term1_efficiency_list, term2_efficiency_list,
                    term3_efficiency_list, gene1_list, gene2_list, gene3_list, prom1_list, prom2_list, term1_list, term2_list, term3_list,
                    rnase1_list, rnase2_list, rnase3_list, i, i_list):

    global gen
    probability = 0
    output_file = "three_genes_replicated.tsv"
    iteration = "best"

    if f_new > f_old:
        #Accepting mutation
        three_genome.recreated_genome(new_pol_strength, new_pol1_strength, new_pol2_strength, new_term1_efficiency,
                                        new_term2_efficiency, new_term3_efficiency, new_gene1_start, new_gene1_stop,
                                        new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop,
                                        promoter1_start, promoter1_stop, promoter2_start, promoter2_stop,
                                        new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop,
                                        new_rnase3_start, new_rnase3_stop, new_term1_start, new_term1_stop,
                                        new_term2_start, new_term2_stop, new_term3_start, new_term3_stop, output_file)
        #Taking in new file and removing unnecessary rows and columns
        nf = pd.read_csv("three_genes_replicated.tsv", header=0, sep='\t')
        nf = edit_new_file(nf)
        sos = sum_of_squares(df, nf)
        all_sos_list.append(sos)
        #Accepts mutation only if sum of squares value decreases
        if sos <= sos_list[-1]:
            output = "best_replicated_gen_" + str(i) + ".tsv"
            poly_list.append(new_pol_strength)
            poly1_list.append(new_pol1_strength)
            poly2_list.append(new_pol2_strength)
            term1_efficiency_list.append(new_term1_efficiency)
            term2_efficiency_list.append(new_term2_efficiency)
            term3_efficiency_list.append(new_term3_efficiency)
            sos_list.append(sos)
            i_list.append(i)
            three_genome.best_recreated_genome(new_pol_strength, new_pol1_strength, new_pol2_strength, new_term1_efficiency,
                                                new_term2_efficiency, new_term3_efficiency, new_gene1_start, new_gene1_stop,
                                                new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop,
                                                promoter1_start, promoter1_stop, promoter2_start, promoter2_stop,
                                                new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop,
                                                new_rnase3_start, new_rnase3_stop, new_term1_start, new_term1_stop,
                                                new_term2_start, new_term2_stop, new_term3_start, new_term3_stop, output)
            #Gets genome coordinates
            gene1_list[0] = new_gene1_start
            gene1_list[1] = new_gene1_stop
            gene2_list[0] = new_gene2_start
            gene2_list[1] = new_gene2_stop
            gene3_list[0] = new_gene3_start
            gene3_list[1] = new_gene3_stop
            prom1_list[0] = promoter1_start
            prom1_list[1] = promoter1_stop
            prom2_list[0] = promoter2_start
            prom2_list[1] = promoter2_stop
            term1_list[0] = new_term1_start
            term1_list[1] = new_term1_stop
            term2_list[0] = new_term2_start
            term2_list[1] = new_term2_stop
            term3_list[0] = new_term3_start
            term3_list[1] = new_term3_stop
            rnase1_list[0] = new_rnase1_start
            rnase1_list[1] = new_rnase1_stop
            rnase2_list[0] = new_rnase2_start
            rnase2_list[1] = new_rnase2_stop
            rnase3_list[0] = new_rnase3_start
            rnase3_list[1] = new_rnase3_stop
            get_coordinates(gene1_list, gene2_list, gene3_list, prom1_list, prom2_list, term1_list, term2_list, term3_list,
                                rnase1_list, rnase2_list, rnase3_list, poly_list, poly1_list, poly2_list, term1_efficiency_list,
                                term2_efficiency_list, term3_efficiency_list, iteration)
            gen+=1
    else:
        #Calculate fitness of new mutation
        probability = calc_fitness(f_new, f_old, Ne)
        f_old = probability
        if probability > random.random():
            #Accepting mutation
            three_genome.recreated_genome(new_pol_strength, new_pol1_strength, new_pol2_strength, new_term1_efficiency,
                                            new_term2_efficiency, new_term3_efficiency, new_gene1_start, new_gene1_stop,
                                            new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop,
                                            promoter1_start, promoter1_stop, promoter2_start, promoter2_stop,
                                            new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop,
                                            new_rnase3_start, new_rnase3_stop, new_term1_start, new_term1_stop,
                                            new_term2_start, new_term2_stop, new_term3_start, new_term3_stop, output_file)
            #Taking in new file and removing unnecessary rows and columns
            nf = pd.read_csv("three_genes_replicated.tsv", header=0, sep='\t')
            nf = edit_new_file(nf)
            sos = sum_of_squares(df, nf)
            all_sos_list.append(sos)
            #Accepts mutation only if sum of squares value decreases
            if sos <= sos_list[-1]:
                output = "best_replicated_gen_" + str(i) + ".tsv"
                poly_list.append(new_pol_strength)
                poly1_list.append(new_pol1_strength)
                poly2_list.append(new_pol2_strength)
                term1_efficiency_list.append(new_term1_efficiency)
                term2_efficiency_list.append(new_term2_efficiency)
                term3_efficiency_list.append(new_term3_efficiency)
                sos_list.append(sos)
                i_list.append(i)
                three_genome.best_recreated_genome(new_pol_strength, new_pol1_strength, new_pol2_strength, new_term1_efficiency,
                                                    new_term2_efficiency, new_term3_efficiency, new_gene1_start, new_gene1_stop,
                                                    new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop,
                                                    promoter1_start, promoter1_stop, promoter2_start, promoter2_stop,
                                                    new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop,
                                                    new_rnase3_start, new_rnase3_stop, new_term1_start, new_term1_stop,
                                                    new_term2_start, new_term2_stop, new_term3_start, new_term3_stop, output)
                #Gets genome coordinates
                gene1_list[0] = new_gene1_start
                gene1_list[1] = new_gene1_stop
                gene2_list[0] = new_gene2_start
                gene2_list[1] = new_gene2_stop
                gene3_list[0] = new_gene3_start
                gene3_list[1] = new_gene3_stop
                prom1_list[0] = promoter1_start
                prom1_list[1] = promoter1_stop
                prom2_list[0] = promoter2_start
                prom2_list[1] = promoter2_stop
                term1_list[0] = new_term1_start
                term1_list[1] = new_term1_stop
                term2_list[0] = new_term2_start
                term2_list[1] = new_term2_stop
                term3_list[0] = new_term3_start
                term3_list[1] = new_term3_stop
                rnase1_list[0] = new_rnase1_start
                rnase1_list[1] = new_rnase1_stop
                rnase2_list[0] = new_rnase2_start
                rnase2_list[1] = new_rnase2_stop
                rnase3_list[0] = new_rnase3_start
                rnase3_list[1] = new_rnase3_stop
                get_coordinates(gene1_list, gene2_list, gene3_list, prom1_list, prom2_list, term1_list, term2_list, term3_list,
                                    rnase1_list, rnase2_list, rnase3_list, poly_list, poly1_list, poly2_list, term1_efficiency_list,
                                    term2_efficiency_list, term3_efficiency_list, iteration)
                gen+=1
    if probability > 0:
        return probability
    else:
        return 0

def alter_poly_strength(poly_list, poly1_list, poly2_list, promoter):

    poly_sigma = 1e5

    if promoter == "promoter":
        #Determining polymerase strength
        poly_eps = np.random.normal(mu, poly_sigma)
        #Determines new polymerase strength to reduce sum of squares value
        pol_strength = poly_list[-1] + poly_eps
        while pol_strength < 0:
            poly_eps = np.random.normal(mu, poly_sigma)
            pol_strength = poly_list[-1] + poly_eps
            if pol_strength > 3e10:
                poly_eps = np.random.normal(mu, poly_sigma)
                pol_strength = poly_list[-1] + poly_eps
        return pol_strength

    if promoter == "promoter1":
        #Determining polymerase strength
        poly_eps = np.random.normal(mu, poly_sigma)
        #Determines new polymerase strength to reduce sum of squares value
        pol1_strength = poly1_list[-1] + poly_eps
        while pol1_strength < 0:
            poly_eps = np.random.normal(mu, poly_sigma)
            pol1_strength = poly1_list[-1] + poly_eps
            if pol1_strength > 3e10:
                poly_eps = np.random.normal(mu, poly_sigma)
                pol1_strength = poly1_list[-1] + poly_eps
        return pol1_strength

    if promoter == "promoter2":
        #Determining polymerase strength
        poly_eps = np.random.normal(mu, poly_sigma)
        #Determines new polymerase strength to reduce sum of squares value
        pol2_strength = poly2_list[-1] + poly_eps
        while pol2_strength < 0:
            poly_eps = np.random.normal(mu, poly_sigma)
            pol2_strength = poly2_list[-1] + poly_eps
            if pol2_strength > 3e10:
                poly_eps = np.random.normal(mu, poly_sigma)
                pol2_strength = poly2_list[-1] + poly_eps
        return pol2_strength

def alter_term_efficiency(term1_efficiency_list, term2_efficiency_list, term3_efficiency_list, terminator):

    term_sigma = 1.0

    if terminator == "terminator1":
        #Determining terminator polymerase efficiency rate
        term_eps = np.random.normal(mu, term_sigma)
        #Determines new terminator efficiency value to reduce sum of squares value
        term1_efficiency = term1_efficiency_list[-1] + term_eps
        while term1_efficiency < 0 or term1_efficiency > 1:
            term_eps = np.random.normal(mu, term_sigma)
            term1_efficiency = term1_efficiency_list[-1] + term_eps
        return term1_efficiency

    if terminator == "terminator2":
        #Determining terminator polymerase efficiency rate
        term_eps = np.random.normal(mu, term_sigma)
        #Determines new terminator efficiency value to reduce sum of squares value
        term2_efficiency = term2_efficiency_list[-1] + term_eps
        while term2_efficiency < 0 or term2_efficiency > 1:
            term_eps = np.random.normal(mu, term_sigma)
            term2_efficiency = term2_efficiency_list[-1] + term_eps
        return term2_efficiency

    if terminator == "terminator3":
        #Determining terminator polymerase efficiency rate
        term_eps = np.random.normal(mu, term_sigma)
        #Determines new terminator efficiency value to reduce sum of squares value
        term3_efficiency = term3_efficiency_list[-1] + term_eps
        while term3_efficiency < 0 or term3_efficiency > 1:
            term_eps = np.random.normal(mu, term_sigma)
            term3_efficiency = term3_efficiency_list[-1] + term_eps
        return term3_efficiency

def add_promoter(prom1_start, prom1_stop, prom2_start, prom2_stop, prom1_start_list, prom1_stop_list, prom2_start_list,
                    prom2_stop_list, term1_start_list, term1_stop_list, term2_start_list, term2_stop_list, rnase1_start_list,
                    rnase1_stop_list, rnase2_start_list, rnase2_stop_list, poly_list, poly1_list, poly2_list):

    region1a_start = 122
    region1a_stop = 132
    region1b_start = 133
    region1b_stop = 143
    region2a_start = 281
    region2a_stop = 291
    region2b_start = 292
    region2b_stop = 302
    pol1_strength = poly1_list[-1]
    pol2_strength = poly2_list[-1]
    promoter_possibilities = ["promoter1", "promoter2"]
    chosen_promoter = random.choice(promoter_possibilities)
    promoter1_slots = ['A', 'B']
    promoter2_slots = ['A', 'B']

    if chosen_promoter == "promoter1":
        promoter = "promoter1"
        #Adding in a promoter between genes 1 and 2
        if (rnase1_start_list[-1] >= region1a_start) and (rnase1_start_list[-1] <= region1a_stop):
            promoter1_slots.remove('A')
        if (term1_start_list[-1] >= region1a_start) and (term1_start_list[-1] <= region1a_stop):
            promoter1_slots.remove('A')
        if (rnase1_start_list[-1] >= region1b_start) and (rnase1_start_list[-1] <= region1b_stop):
            promoter1_slots.remove('B')
        if (term1_start_list[-1] >= region1b_start) and (term1_start_list[-1] <= region1b_stop):
            promoter1_slots.remove('B')
        if promoter1_slots == []:
            prom1_start = prom1_start_list[-1]
            prom1_stop = prom1_stop_list[-1]
        else:
            available_slot = random.choice(promoter1_slots)
            if available_slot == 'A':
                prom1_start = random.randint(region1a_start, 123)
                prom1_stop = prom1_start + 9
            if available_slot == 'B':
                prom1_start = random.randint(region1b_start, 134)
                prom1_stop = prom1_start + 9
        prom1_start_list.append(prom1_start)
        prom1_stop_list.append(prom1_stop)

        pol1_strength = alter_poly_strength(poly_list, poly1_list, poly2_list, promoter)

    if chosen_promoter == "promoter2":
        promoter = "promoter2"
        #Adding promoter between genes 2 and 3
        if (rnase2_start_list[-1] >= region2a_start) and (rnase2_start_list[-1] <= region2a_stop):
            promoter2_slots.remove('A')
        if (term2_start_list[-1] >= region2a_start) and (term2_start_list[-1] <= region2a_stop):
            promoter2_slots.remove('A')
        if (rnase2_start_list[-1] >= region2b_start) and (rnase2_start_list[-1] <= region2b_stop):
            promoter2_slots.remove('B')
        if (term2_start_list[-1] >= region2b_start) and (term2_start_list[-1] <= region2b_stop):
            promoter2_slots.remove('B')
        if promoter2_slots == []:
            prom2_start = prom2_start_list[-1]
            prom2_stop = prom2_stop_list[-1]
        else:
            available_slot = random.choice(promoter2_slots)
            if available_slot == 'A':
                prom2_start = random.randint(region2a_start, 282)
                prom2_stop = prom2_start + 9
            if available_slot == 'B':
                prom2_start = random.randint(region2b_start, 293)
                prom2_stop = prom2_start + 9
        prom2_start_list.append(prom2_start)
        prom2_stop_list.append(prom2_stop)

        pol2_strength = alter_poly_strength(poly_list, poly1_list, poly2_list, promoter)

    prom_list = [pol1_strength, pol2_strength]
    return prom_list

def remove_promoter(promoter1_start, promoter1_stop, promoter2_start, promoter2_stop, prom1_start_list, prom1_stop_list,
                    prom2_start_list, prom2_stop_list):

    promoter_possibilities = ["promoter1", "promoter2"]
    chosen_promoter = random.choice(promoter_possibilities)

    if chosen_promoter == "promoter1":
        #Removing promoter between genes 1 and 2
        promoter1_start = 0
        promoter1_stop = 0
        prom1_start_list.append(promoter1_start)
        prom1_stop_list.append(promoter1_stop)

    if chosen_promoter == "promoter2":
        #Removing promoter between genes 2 and 3
        promoter2_start = 0
        promoter2_stop = 0
        prom2_start_list.append(promoter2_start)
        prom2_stop_list.append(promoter2_stop)

def add_rnase(rnase1_start, rnase1_stop, rnase2_start, rnase2_stop, rnase3_start, rnase3_stop, prom1_start_list,
                prom1_stop_list, prom2_start_list, prom2_stop_list, term1_start_list, term1_stop_list, term2_start_list,
                term2_stop_list, rnase1_start_list, rnase1_stop_list, rnase2_start_list, rnase2_stop_list, rnase3_start_list,
                rnase3_stop_list):

    region1a_start = 122
    region1a_stop = 132
    region1b_start = 133
    region1b_stop = 143
    region1c_start = 144
    region1c_stop = 159
    region2a_start = 281
    region2a_stop = 291
    region2b_start = 292
    region2b_stop = 302
    region2c_start = 303
    region2c_stop = 318
    region_start = 11
    region_stop = 25
    rnase_possibilities = ["rnase1", "rnase2", "rnase3"]
    chosen_rnase = random.choice(rnase_possibilities)
    rnase_slots = ['A', 'B', 'C']

    if chosen_rnase == "rnase1":
        #Adds rnase after first promoter
        rnase1_start = random.randint(region_start, 15)
        rnase1_stop = rnase1_start + 10
        rnase1_start_list.append(rnase1_start)
        rnase1_stop_list.append(rnase1_stop)

    if chosen_rnase == "rnase2":
        #Adds rnase after second promoter
        if (prom1_start_list[-1] >= region1a_start) and (prom1_start_list[-1] <= region1a_stop):
            rnase_slots.remove('A')
        if (term1_start_list[-1] >= region1a_start) and (term1_start_list[-1] <= region1a_stop):
            rnase_slots.remove('A')
        if (prom1_start_list[-1] >= region1b_start) and (prom1_start_list[-1] <= region1b_stop):
            rnase_slots.remove('B')
        if (term1_start_list[-1] >= region1b_start) and (term1_start_list[-1] <= region1b_stop):
            rnase_slots.remove('B')
        if (term1_start_list[-1] >= region1c_start) and (term1_start_list[-1] <= region1c_start):
            rnase_slots.remove('C')
        if rnase_slots == []:
            rnase2_start = rnase2_start_list[-1]
            rnase2_stop = rnase2_stop_list[-1]
        else:
            available_slot = random.choice(rnase_slots)
            if available_slot == 'A':
                rnase2_start = random.randint(region1a_start, 122)
                rnase2_stop = rnase2_start + 10
            if available_slot == 'B':
                rnase2_start = random.randint(region1b_start, 133)
                rnase2_stop = rnase2_start + 10
            if available_slot == 'C':
                rnase2_start = random.randint(region1c_start, 149)
                rnase2_stop = rnase2_start + 10
        rnase2_start_list.append(rnase2_start)
        rnase2_stop_list.append(rnase2_stop)


    if chosen_rnase == "rnase3":
        #Adds rnase after third promoter
        if (prom2_start_list[-1] >= region2a_start) and (prom2_start_list[-1] <= region2a_stop):
            rnase_slots.remove('A')
        if (term2_start_list[-1] >= region2a_start) and (term2_start_list[-1] <= region2a_stop):
            rnase_slots.remove('A')
        if (prom2_start_list[-1] >= region2b_start) and (prom2_start_list[-1] <= region2b_stop):
            rnase_slots.remove('B')
        if (term2_start_list[-1] >= region2b_start) and (term2_start_list[-1] <= region2b_stop):
            rnase_slots.remove('B')
        if (term2_start_list[-1] >= region2c_start) and (term2_start_list[-1] <= region2c_start):
            rnase_slots.remove('C')
        if rnase_slots == []:
            rnase3_start = rnase3_start_list[-1]
            rnase3_stop = rnase3_stop_list[-1]
        else:
            available_slot = random.choice(rnase_slots)
            if available_slot == 'A':
                rnase3_start = random.randint(region2a_start, 281)
                rnase3_stop = rnase3_start + 10
            if available_slot == 'B':
                rnase3_start = random.randint(region2b_start, 292)
                rnase3_stop = rnase3_start + 10
            if available_slot == 'C':
                rnase3_start = random.randint(region2c_start, 308)
                rnase3_stop = rnase3_start + 10
        rnase3_start_list.append(rnase3_start)
        rnase3_stop_list.append(rnase3_stop)

def remove_rnase(rnase1_start, rnase1_stop, rnase2_start, rnase2_stop, rnase3_start, rnase3_stop, rnase1_start_list,
                    rnase1_stop_list, rnase2_start_list, rnase2_stop_list, rnase3_start_list, rnase3_stop_list):

    rnase_possibilities = ["rnase1", "rnase2", "rnase3"]
    chosen_rnase = random.choice(rnase_possibilities)

    if chosen_rnase == "rnase1":
        #Removes rnase after first promoter
        rnase1_start = 0
        rnase1_stop = 0
        rnase1_start_list.append(rnase1_start)
        rnase1_stop_list.append(rnase1_stop)

    if chosen_rnase == "rnase2":
        #Removes rnase after second promoter
        rnase2_start = 0
        rnase2_stop = 0
        rnase2_start_list.append(rnase2_start)
        rnase2_stop_list.append(rnase2_stop)

    if chosen_rnase == "rnase3":
        #Removes rnase after third promoter
        rnase3_start = 0
        rnase3_stop = 0
        rnase3_start_list.append(rnase3_start)
        rnase3_stop_list.append(rnase3_stop)

def add_terminator(term1_start, term1_stop, term2_start, term2_stop, term3_start, term3_stop, prom1_start_list, prom1_stop_list,
                    prom2_start_list, prom2_stop_list, term1_start_list, term1_stop_list, term2_start_list, term2_stop_list,
                    term3_start_list, term3_stop_list, rnase1_start_list, rnase1_stop_list, rnase2_start_list, rnase2_stop_list,
                    term1_efficiency_list, term2_efficiency_list, term3_efficiency_list):

    region1a_start = 122
    region1a_stop = 132
    region1b_start = 133
    region1b_stop = 143
    region1c_start = 144
    region1c_stop = 159
    region2a_start = 281
    region2a_stop = 291
    region2b_start = 292
    region2b_stop = 302
    region2c_start = 303
    region2c_stop = 318
    term1_efficiency = term1_efficiency_list[-1]
    term2_efficiency = term2_efficiency_list[-1]
    term3_efficiency = term3_efficiency_list[-1]
    terminator_possibilities = ["terminator1", "terminator2", "terminator3"]
    chosen_terminator = random.choice(terminator_possibilities)
    terminator_slots = ['A', 'B', 'C']

    if chosen_terminator == "terminator1":
        terminator = "terminator1"
        #Adds terminator after first gene
        if (prom1_start_list[-1] >= region1a_start) and (prom1_start_list[-1] <= region1a_stop):
            terminator_slots.remove('A')
        if (rnase1_start_list[-1] >= region1a_start) and (rnase1_start_list[-1] <= region1a_stop):
            terminator_slots.remove('A')
        if (prom1_start_list[-1] >= region1b_start) and (prom1_start_list[-1] <= region1b_stop):
            terminator_slots.remove('B')
        if (rnase1_start_list[-1] >= region1b_start) and (rnase1_start_list[-1] <= region1b_stop):
            terminator_slots.remove('B')
        if (rnase1_start_list[-1] >= region1c_start) and (rnase1_start_list[-1] <= region1c_start):
            terminator_slots.remove('C')
        if terminator_slots == []:
            term1_start = term1_start_list[-1]
            term1_stop = term1_stop_list[-1]
        else:
            available_slot = random.choice(terminator_slots)
            if available_slot == 'A':
                term1_start = random.randint(region1a_start, 131)
                term1_stop = term1_start + 1
            if available_slot == 'B':
                term1_start = random.randint(region1b_start, 142)
                term1_stop = term1_start + 1
            if available_slot == 'C':
                term1_start = random.randint(region1c_start, 158)
                term1_stop = term1_start + 1
        term1_start_list.append(term1_start)
        term1_stop_list.append(term1_stop)

        term1_efficiency = alter_term_efficiency(term1_efficiency_list, term2_efficiency_list, term3_efficiency_list, terminator)


    if chosen_terminator == "terminator2":
        terminator = "terminator2"
        #Adds terminator after second gene
        if (prom2_start_list[-1] >= region2a_start) and (prom2_start_list[-1] <= region2a_stop):
            terminator_slots.remove('A')
        if (rnase2_start_list[-1] >= region2a_start) and (rnase2_start_list[-1] <= region2a_stop):
            terminator_slots.remove('A')
        if (prom2_start_list[-1] >= region2b_start) and (prom2_start_list[-1] <= region2b_stop):
            terminator_slots.remove('B')
        if (rnase2_start_list[-1] >= region2b_start) and (rnase2_start_list[-1] <= region2b_stop):
            terminator_slots.remove('B')
        if (rnase2_start_list[-1] >= region2c_start) and (rnase2_start_list[-1] <= region2c_start):
            terminator_slots.remove('C')
        if terminator_slots == []:
            term2_start = term2_start_list[-1]
            term2_stop = term2_stop_list[-1]
        else:
            available_slot = random.choice(terminator_slots)
            if available_slot == 'A':
                term2_start = random.randint(region2a_start, 290)
                term2_stop = term2_start + 1
            if available_slot == 'B':
                term2_start = random.randint(region2b_start, 301)
                term2_stop = term2_start + 1
            if available_slot == 'C':
                term2_start = random.randint(region2c_start, 317)
                term2_stop = term2_start + 1
        term2_start_list.append(term2_start)
        term2_stop_list.append(term2_stop)

        term2_efficiency = alter_term_efficiency(term1_efficiency_list, term2_efficiency_list, term3_efficiency_list, terminator)

    if chosen_terminator == "terminator3":
        terminator = "terminator3"
        #Adds terminator after third gene
        term3_start = 449
        term3_stop = 450
        term3_start_list.append(term3_start)
        term3_stop_list.append(term3_stop)

        term3_efficiency = alter_term_efficiency(term1_efficiency_list, term2_efficiency_list, term3_efficiency_list, terminator)

    term_list = [term1_efficiency, term2_efficiency, term3_efficiency]
    return term_list

def remove_terminator(term1_start, term1_stop, term2_start, term2_stop, term3_start, term3_stop, term1_start_list,
                        term1_stop_list, term2_start_list, term2_stop_list, term3_start_list, term3_stop_list):

    terminator_possibilities = ["terminator1", "terminator2", "terminator3"]
    chosen_terminator = random.choice(terminator_possibilities)

    if chosen_terminator == "terminator1":
        #Removes terminator after first gene
        term1_start = 0
        term1_stop = 0
        term1_start_list.append(term1_start)
        term1_stop_list.append(term1_stop)

    if chosen_terminator == "terminator2":
        #Removes terminator after second gene
        term2_start = 0
        term2_stop = 0
        term2_start_list.append(term2_start)
        term2_stop_list.append(term2_stop)

    if chosen_terminator == "terminator3":
        #Removes terminator after third gene
        term3_start = 0
        term3_stop = 0
        term3_start_list.append(term3_start)
        term3_stop_list.append(term3_stop)

def get_coordinates(gene1_list, gene2_list, gene3_list, prom1_list, prom2_list, term1_list, term2_list, term3_list, rnase1_list, rnase2_list,
                        rnase3_list, poly_list, poly1_list, poly2_list, term1_eff_list, term2_eff_list, term3_eff_list, generation):

    best_poly_list = []
    best_poly1_list = []
    best_poly2_list = []
    best_term1_eff_list = []
    best_term2_eff_list = []
    best_term3_eff_list = []

    #Gets coordinates of genome and outputs as tsv file
    if generation == "best":
        best_poly_list.append(poly_list[-1])
        best_poly1_list.append(poly1_list[-1])
        best_poly2_list.append(poly2_list[-1])
        best_term1_eff_list.append(term1_eff_list[-1])
        best_term2_eff_list.append(term2_eff_list[-1])
        best_term3_eff_list.append(term3_eff_list[-1])
        gene1_dataframe = pd.DataFrame(gene1_list, columns=["Gene_Location"])
        export_csv = gene1_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/best_genome/gene1_data.tsv", index=False)
        gene2_dataframe = pd.DataFrame(gene2_list, columns=["Gene_Location"])
        export_csv = gene2_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/best_genome/gene2_data.tsv", index=False)
        gene3_dataframe = pd.DataFrame(gene3_list, columns=["Gene_Location"])
        export_csv = gene3_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/best_genome/gene3_data.tsv", index=False)
        prom1_dataframe = pd.DataFrame(prom1_list, columns=["Promoter_Location"])
        export_csv = prom1_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/best_genome/promoter1_data.tsv", index=False)
        prom2_dataframe = pd.DataFrame(prom2_list, columns=["Promoter_Location"])
        export_csv = prom2_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/best_genome/promoter2_data.tsv", index=False)
        term1_dataframe = pd.DataFrame(term1_list, columns=["Terminator_Location"])
        export_csv = term1_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/best_genome/term1_data.tsv", index=False)
        term2_dataframe = pd.DataFrame(term2_list, columns=["Terminator_Location"])
        export_csv = term2_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/best_genome/term2_data.tsv", index=False)
        term3_dataframe = pd.DataFrame(term3_list, columns=["Terminator_Location"])
        export_csv = term3_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/best_genome/term3_data.tsv", index=False)
        rnase1_dataframe = pd.DataFrame(rnase1_list, columns=["Rnase_Location"])
        export_csv = rnase1_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/best_genome/rnase1_data.tsv", index=False)
        rnase2_dataframe = pd.DataFrame(rnase2_list, columns=["Rnase_Location"])
        export_csv = rnase2_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/best_genome/rnase2_data.tsv", index=False)
        rnase3_dataframe = pd.DataFrame(rnase3_list, columns=["Rnase_Location"])
        export_csv = rnase3_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/best_genome/rnase3_data.tsv", index=False)
        poly_dataframe = pd.DataFrame(best_poly_list, columns=["Polymerase_Rate"])
        export_csv = poly_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/best_genome/poly_data.tsv", index=False)
        poly1_dataframe = pd.DataFrame(best_poly1_list, columns=["Polymerase_Rate"])
        export_csv = poly1_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/best_genome/poly1_data.tsv", index=False)
        poly2_dataframe = pd.DataFrame(best_poly2_list, columns=["Polymerase_Rate"])
        export_csv = poly2_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/best_genome/poly2_data.tsv", index=False)
        term1_efficiency_dataframe = pd.DataFrame(best_term1_eff_list, columns=["Terminator1_Efficiency_Rate"])
        export_csv = term1_efficiency_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/best_genome/term1_efficiency_data.tsv", index=False)
        term2_efficiency_dataframe = pd.DataFrame(best_term2_eff_list, columns=["Terminator2_Efficiency_Rate"])
        export_csv = term2_efficiency_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/best_genome/term2_efficiency_data.tsv", index=False)
        term3_efficiency_dataframe = pd.DataFrame(best_term3_eff_list, columns=["Terminator3_Efficiency_Rate"])
        export_csv = term3_efficiency_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/best_genome/term3_efficiency_data.tsv", index=False)

    if generation == "0":
        best_poly_list.append(poly_list[-1])
        best_poly1_list.append(poly1_list[-1])
        best_poly2_list.append(poly2_list[-1])
        best_term1_eff_list.append(term1_eff_list[-1])
        best_term2_eff_list.append(term2_eff_list[-1])
        best_term3_eff_list.append(term3_eff_list[-1])
        gene1_dataframe = pd.DataFrame(gene1_list, columns=["Gene_Location"])
        export_csv = gene1_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen0/gene1_data.tsv", index=False)
        gene2_dataframe = pd.DataFrame(gene2_list, columns=["Gene_Location"])
        export_csv = gene2_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen0/gene2_data.tsv", index=False)
        gene3_dataframe = pd.DataFrame(gene3_list, columns=["Gene_Location"])
        export_csv = gene3_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen0/gene3_data.tsv", index=False)
        prom1_dataframe = pd.DataFrame(prom1_list, columns=["Promoter_Location"])
        export_csv = prom1_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen0/promoter1_data.tsv", index=False)
        prom2_dataframe = pd.DataFrame(prom2_list, columns=["Promoter_Location"])
        export_csv = prom2_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen0/promoter2_data.tsv", index=False)
        term1_dataframe = pd.DataFrame(term1_list, columns=["Terminator_Location"])
        export_csv = term1_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen0/term1_data.tsv", index=False)
        term2_dataframe = pd.DataFrame(term2_list, columns=["Terminator_Location"])
        export_csv = term2_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen0/term2_data.tsv", index=False)
        term3_dataframe = pd.DataFrame(term3_list, columns=["Terminator_Location"])
        export_csv = term3_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen0/term3_data.tsv", index=False)
        rnase1_dataframe = pd.DataFrame(rnase1_list, columns=["Rnase_Location"])
        export_csv = rnase1_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen0/rnase1_data.tsv", index=False)
        rnase2_dataframe = pd.DataFrame(rnase2_list, columns=["Rnase_Location"])
        export_csv = rnase2_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen0/rnase2_data.tsv", index=False)
        rnase3_dataframe = pd.DataFrame(rnase3_list, columns=["Rnase_Location"])
        export_csv = rnase3_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen0/rnase3_data.tsv", index=False)
        poly_dataframe = pd.DataFrame(best_poly_list, columns=["Polymerase_Rate"])
        export_csv = poly_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen0/poly_data.tsv", index=False)
        poly1_dataframe = pd.DataFrame(best_poly1_list, columns=["Polymerase_Rate"])
        export_csv = poly1_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen0/poly1_data.tsv", index=False)
        poly2_dataframe = pd.DataFrame(best_poly2_list, columns=["Polymerase_Rate"])
        export_csv = poly2_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen0/poly2_data.tsv", index=False)
        term1_efficiency_dataframe = pd.DataFrame(best_term1_eff_list, columns=["Terminator1_Efficiency_Rate"])
        export_csv = term1_efficiency_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen0/term1_efficiency_data.tsv", index=False)
        term2_efficiency_dataframe = pd.DataFrame(best_term2_eff_list, columns=["Terminator2_Efficiency_Rate"])
        export_csv = term2_efficiency_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen0/term2_efficiency_data.tsv", index=False)
        term3_efficiency_dataframe = pd.DataFrame(best_term3_eff_list, columns=["Terminator3_Efficiency_Rate"])
        export_csv = term3_efficiency_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen0/term3_efficiency_data.tsv", index=False)

    if generation == "5000":
        best_poly_list.append(poly_list[-1])
        best_poly1_list.append(poly1_list[-1])
        best_poly2_list.append(poly2_list[-1])
        best_term1_eff_list.append(term1_eff_list[-1])
        best_term2_eff_list.append(term2_eff_list[-1])
        best_term3_eff_list.append(term3_eff_list[-1])
        gene1_dataframe = pd.DataFrame(gene1_list, columns=["Gene_Location"])
        export_csv = gene1_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen5000/gene1_data.tsv", index=False)
        gene2_dataframe = pd.DataFrame(gene2_list, columns=["Gene_Location"])
        export_csv = gene2_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen5000/gene2_data.tsv", index=False)
        gene3_dataframe = pd.DataFrame(gene3_list, columns=["Gene_Location"])
        export_csv = gene3_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen5000/gene3_data.tsv", index=False)
        prom1_dataframe = pd.DataFrame(prom1_list, columns=["Promoter_Location"])
        export_csv = prom1_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen5000/promoter1_data.tsv", index=False)
        prom2_dataframe = pd.DataFrame(prom2_list, columns=["Promoter_Location"])
        export_csv = prom2_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen5000/promoter2_data.tsv", index=False)
        term1_dataframe = pd.DataFrame(term1_list, columns=["Terminator_Location"])
        export_csv = term1_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen5000/term1_data.tsv", index=False)
        term2_dataframe = pd.DataFrame(term2_list, columns=["Terminator_Location"])
        export_csv = term2_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen5000/term2_data.tsv", index=False)
        term3_dataframe = pd.DataFrame(term3_list, columns=["Terminator_Location"])
        export_csv = term3_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen5000/term3_data.tsv", index=False)
        rnase1_dataframe = pd.DataFrame(rnase1_list, columns=["Rnase_Location"])
        export_csv = rnase1_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen5000/rnase1_data.tsv", index=False)
        rnase2_dataframe = pd.DataFrame(rnase2_list, columns=["Rnase_Location"])
        export_csv = rnase2_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen5000/rnase2_data.tsv", index=False)
        rnase3_dataframe = pd.DataFrame(rnase3_list, columns=["Rnase_Location"])
        export_csv = rnase3_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen5000/rnase3_data.tsv", index=False)
        poly_dataframe = pd.DataFrame(best_poly_list, columns=["Polymerase_Rate"])
        export_csv = poly_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen5000/poly_data.tsv", index=False)
        poly1_dataframe = pd.DataFrame(best_poly1_list, columns=["Polymerase_Rate"])
        export_csv = poly1_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen5000/poly1_data.tsv", index=False)
        poly2_dataframe = pd.DataFrame(best_poly2_list, columns=["Polymerase_Rate"])
        export_csv = poly2_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen5000/poly2_data.tsv", index=False)
        term1_efficiency_dataframe = pd.DataFrame(best_term1_eff_list, columns=["Terminator1_Efficiency_Rate"])
        export_csv = term1_efficiency_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen5000/term1_efficiency_data.tsv", index=False)
        term2_efficiency_dataframe = pd.DataFrame(best_term2_eff_list, columns=["Terminator2_Efficiency_Rate"])
        export_csv = term2_efficiency_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen5000/term2_efficiency_data.tsv", index=False)
        term3_efficiency_dataframe = pd.DataFrame(best_term3_eff_list, columns=["Terminator3_Efficiency_Rate"])
        export_csv = term3_efficiency_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen5000/term3_efficiency_data.tsv", index=False)

    if generation == "10000":
        best_poly_list.append(poly_list[-1])
        best_poly1_list.append(poly1_list[-1])
        best_poly2_list.append(poly2_list[-1])
        best_term1_eff_list.append(term1_eff_list[-1])
        best_term2_eff_list.append(term2_eff_list[-1])
        best_term3_eff_list.append(term3_eff_list[-1])
        gene1_dataframe = pd.DataFrame(gene1_list, columns=["Gene_Location"])
        export_csv = gene1_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen10000/gene1_data.tsv", index=False)
        gene2_dataframe = pd.DataFrame(gene2_list, columns=["Gene_Location"])
        export_csv = gene2_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen10000/gene2_data.tsv", index=False)
        gene3_dataframe = pd.DataFrame(gene3_list, columns=["Gene_Location"])
        export_csv = gene3_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen10000/gene3_data.tsv", index=False)
        prom1_dataframe = pd.DataFrame(prom1_list, columns=["Promoter_Location"])
        export_csv = prom1_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen10000/promoter1_data.tsv", index=False)
        prom2_dataframe = pd.DataFrame(prom2_list, columns=["Promoter_Location"])
        export_csv = prom2_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen10000/promoter2_data.tsv", index=False)
        term1_dataframe = pd.DataFrame(term1_list, columns=["Terminator_Location"])
        export_csv = term1_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen10000/term1_data.tsv", index=False)
        term2_dataframe = pd.DataFrame(term2_list, columns=["Terminator_Location"])
        export_csv = term2_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen10000/term2_data.tsv", index=False)
        term3_dataframe = pd.DataFrame(term3_list, columns=["Terminator_Location"])
        export_csv = term3_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen10000/term3_data.tsv", index=False)
        rnase1_dataframe = pd.DataFrame(rnase1_list, columns=["Rnase_Location"])
        export_csv = rnase1_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen10000/rnase1_data.tsv", index=False)
        rnase2_dataframe = pd.DataFrame(rnase2_list, columns=["Rnase_Location"])
        export_csv = rnase2_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen10000/rnase2_data.tsv", index=False)
        rnase3_dataframe = pd.DataFrame(rnase3_list, columns=["Rnase_Location"])
        export_csv = rnase3_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen10000/rnase3_data.tsv", index=False)
        poly_dataframe = pd.DataFrame(best_poly_list, columns=["Polymerase_Rate"])
        export_csv = poly_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen10000/poly_data.tsv", index=False)
        poly1_dataframe = pd.DataFrame(best_poly1_list, columns=["Polymerase_Rate"])
        export_csv = poly1_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen10000/poly1_data.tsv", index=False)
        poly2_dataframe = pd.DataFrame(best_poly2_list, columns=["Polymerase_Rate"])
        export_csv = poly2_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen10000/poly2_data.tsv", index=False)
        term1_efficiency_dataframe = pd.DataFrame(best_term1_eff_list, columns=["Terminator1_Efficiency_Rate"])
        export_csv = term1_efficiency_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen10000/term1_efficiency_data.tsv", index=False)
        term2_efficiency_dataframe = pd.DataFrame(best_term2_eff_list, columns=["Terminator2_Efficiency_Rate"])
        export_csv = term2_efficiency_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen10000/term2_efficiency_data.tsv", index=False)
        term3_efficiency_dataframe = pd.DataFrame(best_term3_eff_list, columns=["Terminator3_Efficiency_Rate"])
        export_csv = term3_efficiency_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen10000/term3_efficiency_data.tsv", index=False)

#Class containing genome for simulation with new polymerase strength
class three_genome:

    def __init__(self):
        self.pol_strength = 4e10

    def recreated_genome(pol_strength, pol1_strength, pol2_strength, term1_efficiency, term2_efficiency,
                         term3_efficiency, gene1_start, gene1_stop, gene2_start, gene2_stop, gene3_start,
                         gene3_stop, prom1_start, prom1_stop, prom2_start, prom2_stop, rnase1_start,
                         rnase1_stop, rnase2_start, rnase2_stop, rnase3_start, rnase3_stop, term1_start,
                         term1_stop, term2_start, term2_stop, term3_start, term3_stop, file_name):

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
                             interactions={"rnapol": pol_strength})
        if prom1_start > 0:
            plasmid.add_promoter(name="p2", start=prom1_start, stop=prom1_stop,
                                 interactions={"rnapol": pol1_strength})
        if prom2_start > 0:
            plasmid.add_promoter(name="p3", start=prom2_start, stop=prom2_stop,
                                 interactions={"rnapol": pol2_strength})
        if rnase1_start > 0:
            plasmid.add_rnase_site(start=rnase1_start, stop=rnase1_stop)
        if rnase2_start > 0:
            plasmid.add_rnase_site(start=rnase2_start, stop=rnase2_stop)
        if rnase3_start > 0:
            plasmid.add_rnase_site(start=rnase3_start, stop=rnase3_stop)
        if term1_start > 0:
            plasmid.add_terminator(name="t1", start=term1_start, stop=term1_stop,
                                   efficiency={"rnapol": term1_efficiency})
        if term2_start > 0:
            plasmid.add_terminator(name="t2", start=term2_start, stop=term2_stop,
                                   efficiency={"rnapol": term2_efficiency})
        if term3_start > 0:
            plasmid.add_terminator(name="t3", start=term3_start, stop=term3_stop,
                                   efficiency={"rnapol": term3_efficiency})
        plasmid.add_gene(name="proteinX", start=gene1_start, stop=gene1_stop,
                         rbs_start=(gene1_start-15), rbs_stop=gene1_start, rbs_strength=1e7)

        plasmid.add_gene(name="proteinY", start=gene2_start, stop=gene2_stop,
                         rbs_start=(gene2_start-15), rbs_stop=gene2_start, rbs_strength=1e7)
        plasmid.add_gene(name="proteinZ", start=gene3_start, stop=gene3_stop,
                         rbs_start=(gene3_start-15), rbs_stop=gene3_start, rbs_strength=1e7)
        sim.register_genome(plasmid)
        sim.simulate(time_limit=240, time_step=1,
                     output = file_name)

    def best_recreated_genome(pol_strength, pol1_strength, pol2_strength, term1_efficiency, term2_efficiency,
                                term3_efficiency, gene1_start, gene1_stop, gene2_start, gene2_stop, gene3_start,
                                gene3_stop, prom1_start, prom1_stop, prom2_start, prom2_stop, rnase1_start,
                                rnase1_stop, rnase2_start, rnase2_stop, rnase3_start, rnase3_stop, term1_start,
                                term1_stop, term2_start, term2_stop, term3_start, term3_stop, file_name):

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
                             interactions={"rnapol": pol_strength})
        if prom1_start > 0:
            plasmid.add_promoter(name="p2", start=prom1_start, stop=prom1_stop,
                                 interactions={"rnapol": pol1_strength})
        if prom2_start > 0:
            plasmid.add_promoter(name="p3", start=prom2_start, stop=prom2_stop,
                                 interactions={"rnapol": pol2_strength})
        if rnase1_start > 0:
            plasmid.add_rnase_site(start=rnase1_start, stop=rnase1_stop)
        if rnase2_start > 0:
            plasmid.add_rnase_site(start=rnase2_start, stop=rnase2_stop)
        if rnase3_start > 0:
            plasmid.add_rnase_site(start=rnase3_start, stop=rnase3_stop)
        if term1_start > 0:
            plasmid.add_terminator(name="t1", start=term1_start, stop=term1_stop,
                                   efficiency={"rnapol": term1_efficiency})
        if term2_start > 0:
            plasmid.add_terminator(name="t2", start=term2_start, stop=term2_stop,
                                   efficiency={"rnapol": term2_efficiency})
        if term3_start > 0:
            plasmid.add_terminator(name="t3", start=term3_start, stop=term3_stop,
                                   efficiency={"rnapol": term3_efficiency})
        plasmid.add_gene(name="proteinX", start=gene1_start, stop=gene1_stop,
                         rbs_start=(gene1_start-15), rbs_stop=gene1_start, rbs_strength=1e7)

        plasmid.add_gene(name="proteinY", start=gene2_start, stop=gene2_stop,
                         rbs_start=(gene2_start-15), rbs_stop=gene2_start, rbs_strength=1e7)
        plasmid.add_gene(name="proteinZ", start=gene3_start, stop=gene3_stop,
                         rbs_start=(gene3_start-15), rbs_stop=gene3_start, rbs_strength=1e7)
        sim.register_genome(plasmid)
        sim.simulate(time_limit=240, time_step=1,
                     output = file_name)


if __name__ == '__main__':
    main()
