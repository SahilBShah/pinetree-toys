import pinetree as pt
import pandas
import random
import numpy as np
import string
import filecmp
import math

#Need to add in the following functions: change numeric values according to biological range, check to see if we need to restrict movement after gene sequence has been rearranged, figure out what to do with the global variables

def main():

    global term_list
    global all_term_list
    global poly_list
    global sos_list
    global all_sos_list
    global all_poly_list
    global Ne
    global gen
    global mu
    global sigma
    global gene1_start_list
    global gene2_start_list
    global gene3_start_list
    global gene1_stop_list
    global gene2_stop_list
    global gene3_stop_list
    global term1_start_list
    global term1_stop_list
    global term2_start_list
    global term2_stop_list
    global term3_start_list
    global term3_stop_list
    global rnase1_start_list
    global rnase1_stop_list
    global rnase2_start_list
    global rnase2_stop_list
    global rnase3_start_list
    global rnase3_stop_list
    global prom1_start_list
    global prom1_stop_list
    global prom2_start_list
    global prom2_stop_list
    global new_gene1_start
    global new_gene1_stop
    global new_gene2_start
    global new_gene2_stop
    global new_gene3_start
    global new_gene3_stop
    global promoter1_start
    global promoter1_stop
    global promoter2_start
    global promoter2_stop
    global new_prom1_start
    global new_prom1_stop
    global new_prom2_start
    global new_prom2_stop
    global new_rnase1_start
    global new_rnase1_stop
    global new_rnase2_start
    global new_rnase2_stop
    global new_rnase3_start
    global new_rnase3_stop
    global new_term1_start
    global new_term1_stop
    global new_term2_start
    global new_term2_stop
    global new_term3_start
    global new_term3_stop
    sos_list = [60000]
    poly_list = []
    all_poly_list = []
    all_sos_list = []
    term_list = []
    all_term_list = []
    gene1_start_list = []
    gene2_start_list = []
    gene3_start_list = []
    gene1_stop_list = []
    gene2_stop_list = []
    gene3_stop_list = []
    term1_start_list = [0]
    term1_stop_list = [0]
    term2_start_list = [0]
    term2_stop_list = [0]
    term3_start_list = [0]
    term3_stop_list = [0]
    rnase1_start_list = [0]
    rnase1_stop_list = [0]
    rnase2_start_list = [0]
    rnase2_stop_list = [0]
    rnase3_start_list = [0]
    rnase3_stop_list = [0]
    prom1_start_list = [1]
    prom1_stop_list = [10]
    prom2_start_list = [1]
    prom2_stop_list = [10]
    f_old = 1.0
    mu = 0.0
    sigma = 1.0
    Ne = 10
    i = 0
    gen = 0
    new_gene1_start = 26
    new_gene1_stop = 148
    new_gene2_start = 178
    new_gene2_stop = 300
    new_gene3_start = 315
    new_gene3_stop = 448
    promoter1_start = 0
    promoter1_stop = 0
    promoter2_start = 0
    promoter2_stop = 0
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
    df = pandas.read_table(input("Enter tsv file name: "), delim_whitespace=True, header=0)
    df = edit_target_file(df)

    #Evolution program
    new_pol_strength = random.randint(1e10, 5e10)
    new_term_efficiency = random.uniform(0.0, 1.0)
    poly_list.append(new_pol_strength)
    all_poly_list.append(new_pol_strength)
    term_list.append(new_term_efficiency)
    all_term_list.append(new_term_efficiency)

    #Determining gene sizes
    new_gene1_start = random.randint(25, 449)
    new_gene1_stop = random.randint(26, new_gene2_start-15)
    gene1_stop_list.append(new_gene1_stop)
    while new_gene1_start > new_gene1_stop:
        new_gene1_start = random.randint(25, gene1_stop_list[-1])
        new_gene1_stop = random.randint(26, new_gene2_start-15)
    new_gene2_start = random.randint(new_gene1_stop+15, 400)
    new_gene2_stop = random.randint(new_gene2_start+1, new_gene3_stop-15)
    gene2_stop_list.append(new_gene2_stop)
    while new_gene2_start > new_gene2_stop:
        new_gene2_start = random.randint(new_gene1_stop+15, gene2_stop_list[-1])
        new_gene2_stop = random.randint(new_gene2_start+1, new_gene3_stop-15)
    new_gene3_start = random.randint(new_gene2_stop+15, 449)
    new_gene3_stop = random.randint(new_gene3_start+1, 450)
    gene3_stop_list.append(new_gene3_stop)
    while new_gene3_start > new_gene3_stop:
        new_gene3_start = random.randint(new_gene2_stop+15, gene3_stop_list[-1])
        new_gene3_stop = random.randint(new_gene3_start+1, 450)
    gene1_start_list.append(new_gene1_start)
    gene2_start_list.append(new_gene2_start)
    gene3_start_list.append(new_gene3_start)


    possibilities = [alter_poly_strength(new_pol_strength),
                        alter_term_efficiency(new_term_efficiency),
                        alter_gene_length(new_gene1_start, new_gene1_stop, new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop),
                        add_promoter(new_prom1_start, new_prom1_stop, new_prom2_start, new_prom2_stop),
                        remove_promoter(promoter1_start, promoter1_stop, promoter2_start, promoter2_stop),
                        add_ranse(promoter1_start, promoter2_start),
                        add_terminator(new_gene2_start, new_gene1_stop, new_gene3_start, new_gene2_stop, new_gene3_stop),
                        remove_terminator(new_term1_start, new_term1_stop, new_term2_start, new_term2_stop, new_term3_start, new_term3_stop),
                        remove_rnase(new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop, new_rnase3_start, new_rnase3_stop)]

    while i < 100000:

        random.choice(possibilities)
        eps = np.random.normal(mu, sigma)
        f_new = f_old * (1.0 + eps)
        test = accept_mutation(df, f_old, f_new, new_pol_strength, new_term_efficiency, new_gene1_start, new_gene1_stop, new_gene2_start,
                                    new_gene2_stop, new_gene3_start, new_gene3_stop, promoter1_start, promoter1_stop, promoter2_start,
                                    promoter2_stop, new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop, new_rnase3_start,
                                    new_rnase3_stop, new_term1_start, new_term1_stop, new_term2_start, new_term2_stop, new_term3_start,
                                    new_term3_stop)
        if test > 0:
            f_old = test
        if sos_list[-1] == 0:
            break

        i+=1
        print("i = ", i)


    #Exported tsv files
    all_sos_dataframe = pandas.DataFrame(all_sos_list, columns=["Sum of Squares"])
    export_csv = all_sos_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/all_sos_data.tsv", index=False)
    sos_dataframe = pandas.DataFrame(sos_list, columns=["Sum of Squares"])
    export_csv = sos_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/sos_data.tsv", index=False)
    all_poly_dataframe = pandas.DataFrame(all_poly_list, columns=["Polymerase Rate"])
    export_csv = all_poly_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/all_poly_data.tsv", index=False)
    poly_dataframe = pandas.DataFrame(poly_list, columns=["Polymerase Rate"])
    export_csv = poly_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/poly_data.tsv", index=False)
    all_term_dataframe = pandas.DataFrame(all_term_list, columns=["Terminator Efficiency Rate"])
    export_csv = all_term_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/all_term_data.tsv", index=False)
    term_dataframe = pandas.DataFrame(term_list, columns=["Terminator Efficiency Rate"])
    export_csv = term_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/term_data.tsv", index=False)
    gene1_start_dataframe = pandas.DataFrame(gene1_start_list, columns=["Gene Start"])
    export_csv = gene1_start_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/gene1_start_data.tsv", index=False)
    gene2_start_dataframe = pandas.DataFrame(gene2_start_list, columns=["Gene Start"])
    export_csv = gene2_start_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/gene2_start_data.tsv", index=False)
    gene3_start_dataframe = pandas.DataFrame(gene3_start_list, columns=["Gene Start"])
    export_csv = gene3_start_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/gene3_start_data.tsv", index=False)
    gene1_stop_dataframe = pandas.DataFrame(gene1_stop_list, columns=["Gene Stop"])
    export_csv = gene1_stop_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/gene1_stop_data.tsv", index=False)
    gene2_stop_dataframe = pandas.DataFrame(gene2_stop_list, columns=["Gene Stop"])
    export_csv = gene2_stop_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/gene2_stop_data.tsv", index=False)
    gene3_stop_dataframe = pandas.DataFrame(gene3_stop_list, columns=["Gene Stop"])
    export_csv = gene3_stop_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/gene3_stop_data.tsv", index=False)
    print("Final Polymerase Rate: ", poly_list[-1])
    print("Final Terminator Efficiency Rate: ", term_list[-1])
    print("Gene 1 starts at: ", new_gene1_start, "and stops at: ", new_gene1_stop)
    print("Gene 1 Ribosome Binding Site starts at: ", new_gene1_start-15, "and stops at: ", new_gene1_start)
    print("Gene 1 terminator starts at: ", new_term1_start, "and stop at: ", new_term1_stop)
    print("Gene 2 starts at: ", new_gene2_start, "and stop at: ", new_gene2_stop)
    print("Gene 2 Ribosome Binding Site starts at: ", new_gene2_start-15, "and stop at: ", new_gene2_start)
    print("Gene 2 terminator starts at: ", new_term2_start, "and stop at: ", new_term2_stop)
    print("Gene 3 starts at: ", new_gene3_start, "and stop at: ", new_gene3_stop)
    print("Gene 3 Ribosome Binding Site starts at: ", new_gene3_start-15, "and stop at: ", new_gene3_start)
    print("Gene 3 terminator starts at: ", new_term3_start, "and stop at: ", new_term3_stop)
    print("Promoter 1 starts at: ", 1, "and stops at: ", 10)
    print("Added Promoter 2 starts at: ", promoter1_start, "and stop at: ", promoter1_stop)
    print("Added Promoter 3 starts at: ", promoter2_start, "and stop at: ", promoter2_stop)
    print("Rnase site 1 starts at: ", new_rnase1_start, "and stop at: ", new_rnase1_stop)
    print("Rnase site 2 starts at: ", new_rnase2_start, "and stop at: ", new_rnase2_stop)
    print("Rnase site 3 starts at: ", new_rnase3_start, "and stop at: ", new_rnase3_stop)
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
        p =((1-pow((xi/xj), 2)) /(1-pow((xi/xj), (2 * float(N)))))
    except OverflowError as e:
        p = 0.0
    return (p)

#Minimizes distance between observed values and regression line from given data
def sum_of_squares(target_file, new_file):

    sos = 0.0
    for row, row2 in zip(new_file.iterrows(), target_file.iterrows()):
        sos = sos + (row[1][1] - row2[1][1]) ** 2
    return sos

#Removes unnecessary rows and columns in produced file
def edit_new_file(new_file):

    new_file = new_file.drop("time", axis=1)
    new_file = new_file.drop(columns="protein", axis=1)
    new_file = new_file.drop(columns="ribo_density", axis=1)
    new_file = new_file[new_file.species != '__proteinX_rbs']
    new_file = new_file[new_file.species != '__proteinY_rbs']
    new_file = new_file[new_file.species != '__proteinZ_rbs']
    new_file = new_file[new_file.species != '__ribosome']
    new_file = new_file[new_file.species != '__rnase_site']
    new_file = new_file[new_file.species != '__rnase_site_ext']
    new_file = new_file[new_file.species != 'rnapol']
    new_file = new_file[new_file.species != 'p1']
    new_file = new_file[new_file.species != 'p2']
    return new_file

#Removes unnecessary rows and columns in target file
def edit_target_file(target_file):

    target_file = target_file.drop("time", axis=1)
    target_file = target_file.drop(columns="protein", axis=1)
    target_file = target_file.drop(columns="ribo_density", axis=1)
    target_file = target_file[target_file.species != '__proteinX_rbs']
    target_file = target_file[target_file.species != '__proteinY_rbs']
    target_file = target_file[target_file.species != '__proteinZ_rbs']
    target_file = target_file[target_file.species != '__ribosome']
    target_file = target_file[target_file.species != '__rnase_site']
    target_file = target_file[target_file.species != '__rnase_site_ext']
    target_file = target_file[target_file.species != 'rnapol']
    target_file = target_file[target_file.species != 'p1']
    target_file = target_file[target_file.species != 'p2']
    return target_file

def accept_mutation(df, f_old, f_new, new_pol_strength, new_term_efficiency, new_gene1_start, new_gene1_stop, new_gene2_start,
                    new_gene2_stop, new_gene3_start, new_gene3_stop, promoter1_start, promoter1_stop, promoter2_start,
                    promoter2_stop, new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop, new_rnase3_start,
                    new_rnase3_stop, new_term1_start, new_term1_stop, new_term2_start, new_term2_stop, new_term3_start,
                    new_term3_stop):

    global gen
    probability = 0

    if f_new > f_old:
        #Accepting mutation
        three_genome.recreated_genome(new_pol_strength, new_term_efficiency, new_gene1_start, new_gene1_stop,
                                        new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop,
                                        promoter1_start, promoter1_stop, promoter2_start, promoter2_stop,
                                        new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop,
                                        new_rnase3_start, new_rnase3_stop, new_term1_start, new_term1_stop,
                                        new_term2_start, new_term2_stop, new_term3_start, new_term3_stop)
        #Taking in new file and removing unnecessary rows and columns
        nf = pandas.read_table("three_genes_replicated.tsv", delim_whitespace=True, header=0)
        nf = edit_new_file(nf)
        sos = sum_of_squares(df, nf)
        all_sos_list.append(sos)
        #Accepts mutation only if sum of squares value decreases
        if sos <= sos_list[-1]:
            poly_list.append(new_pol_strength)
            term_list.append(new_term_efficiency)
            sos_list.append(sos)
            three_genome.best_recreated_genome(new_pol_strength, new_term_efficiency, new_gene1_start, new_gene1_stop,
                                                new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop,
                                                promoter1_start, promoter1_stop, promoter2_start, promoter2_stop,
                                                new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop,
                                                new_rnase3_start, new_rnase3_stop, new_term1_start, new_term1_stop,
                                                new_term2_start, new_term2_stop, new_term3_start, new_term3_stop)
            gen+=1
    else:
        #Calculate fitness of new mutation
        probability = calc_fitness(f_new, f_old, Ne)
        f_old = probability
        if probability > random.random():
            #Accepting mutation
            three_genome.recreated_genome(new_pol_strength, new_term_efficiency, new_gene1_start, new_gene1_stop,
                                            new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop,
                                            promoter1_start, promoter1_stop, promoter2_start, promoter2_stop,
                                            new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop,
                                            new_rnase3_start, new_rnase3_stop, new_term1_start, new_term1_stop,
                                            new_term2_start, new_term2_stop, new_term3_start, new_term3_stop)
            #Taking in new file and removing unnecessary rows and columns
            nf = pandas.read_table("three_genes_replicated.tsv", delim_whitespace=True, header=0)
            nf = edit_new_file(nf)
            sos = sum_of_squares(df, nf)
            all_sos_list.append(sos)
            #Accepts mutation only if sum of squares value decreases
            if sos <= sos_list[-1]:
                poly_list.append(new_pol_strength)
                term_list.append(new_term_efficiency)
                sos_list.append(sos)
                three_genome.best_recreated_genome(new_pol_strength, new_term_efficiency, new_gene1_start, new_gene1_stop,
                                                    new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop,
                                                    promoter1_start, promoter1_stop, promoter2_start, promoter2_stop,
                                                    new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop,
                                                    new_rnase3_start, new_rnase3_stop, new_term1_start, new_term1_stop,
                                                    new_term2_start, new_term2_stop, new_term3_start, new_term3_stop)
                gen+=1
    if probability > 0:
        return probability
    else:
        return 0

def alter_poly_strength(new_pol_strength):

    poly_sigma = 1e10
    #Determining polymerase strength
    poly_eps = np.random.normal(mu, poly_sigma)
    #Determines new polymerase strength to reduce sum of squares value
    new_pol_strength = poly_list[-1] + poly_eps
    while new_pol_strength < 0:
        poly_eps = np.random.normal(mu, poly_sigma)
        new_pol_strength = poly_list[-1] + poly_eps
        if new_pol_strength in all_poly_list:
            poly_eps = np.random.normal(mu, poly_sigma)
            new_pol_strength = poly_list[-1] + poly_eps
    if new_pol_strength in all_poly_list:
        poly_eps = np.random.normal(mu, poly_sigma)
        new_pol_strength = poly_list[-1] + poly_eps
    all_poly_list.append(new_pol_strength)

def alter_term_efficiency(new_term_efficiency):

    term_sigma = 1.0
    #Determining terminator polymerase efficiency rate
    term_eps = np.random.normal(mu, term_sigma)
    #Determines new terminator efficiency value to reduce sum of squares value
    new_term_efficiency = term_list[-1] + term_eps
    while new_term_efficiency < 0 or new_term_efficiency > 1:
        term_eps = np.random.normal(mu, term_sigma)
        new_term_efficiency = term_list[-1] + term_eps
        if new_term_efficiency in all_term_list:
            term_eps = np.random.normal(mu, term_sigma)
            new_term_efficiency = term_list[-1] + term_eps
    if new_term_efficiency in all_term_list:
        term_eps = np.random.normal(mu, term_sigma)
        new_term_efficiency = term_list[-1] + term_eps
    all_term_list.append(new_term_efficiency)

def alter_gene_length(new_gene1_start, new_gene1_stop, new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop):

    gene_sigma = 10.0
    gene_possibilities = ["gene1", "gene2", "gene3"]
    chosen_gene = random.choice(gene_possibilities)
    if chosen_gene == "gene1":
        #Changing gene lengths
        gene1_eps = np.random.normal(mu, gene_sigma)
        #Alter gene sizes
        new_gene1_start += round(gene1_eps)
        new_gene1_stop += round(gene1_eps)

        #If gene size is larger or smaller than genome, redo
        while new_gene3_stop > 450 or new_gene1_start < 26:
            gene1_eps = np.random.normal(mu, gene_sigma)
            new_gene1_start = gene1_start_list[-1] + round(gene1_eps)
            new_gene1_stop = gene1_stop_list[-1] + round(gene1_eps)
        while (new_gene1_start > term1_start_list[-1] and new_gene1_start < term1_stop_list[-1]) or (new_gene1_stop > term1_start_list[-1] and new_gene1_stop < term1_stop_list[-1]):
            gene1_eps = np.random.normal(mu, gene_sigma)
            new_gene1_start = gene1_start_list[-1] + round(gene1_eps)
            new_gene1_stop = gene1_stop_list[-1] + round(gene1_eps)
        while (new_gene1_start > prom1_start_list[-1] and new_gene1_start < prom1_stop_list[-1]) or (new_gene1_stop > prom1_start_list[-1] and new_gene1_stop < prom1_stop_list[-1]):
            gene1_eps = np.random.normal(mu, gene_sigma)
            new_gene1_start = gene1_start_list[-1] + round(gene1_eps)
            new_gene1_stop = gene1_stop_list[-1] + round(gene1_eps)
        while (new_gene1_start > rnase1_start_list[-1] and new_gene1_start < rnase1_stop_list[-1]) or (new_gene1_stop > rnase1_start_list[-1] and new_gene1_stop < rnase1_stop_list[-1]):
            gene1_eps = np.random.normal(mu, gene_sigma)
            new_gene1_start = gene1_start_list[-1] + round(gene1_eps)
            new_gene1_stop = gene1_stop_list[-1] + round(gene1_eps)
        #if (new_gene1_start > prom1_start_list[-1] and new_gene1_start < prom1_stop_list[-1]) or (new_gene1_stop > prom1_start_list[-1] and new_gene1_stop < prom1_stop_list[-1]) or (new_gene1_start > prom1_start_list[-1] and new_gene1_start < prom1_stop_list[-1]) or (new_gene1_stop > prom1_start_list[-1] and new_gene1_stop < prom1_stop_list[-1]) or (new_gene1_start > rnase1_start_list[-1] and new_gene1_start < rnase1_stop_list[-1]) or (new_gene1_stop > rnase1_start_list[-1] and new_gene1_stop < rnase1_stop_list[-1]) or new_gene3_stop > 450 or new_gene1_start < 26:
            #alter_gene1_length(new_gene1_start, new_gene1_stop)

        gene1_start_list.append(new_gene1_start)
        gene1_stop_list.append(new_gene1_stop)

    if chosen_gene == "gene2":
        #Changing gene lengths
        gene2_eps = np.random.normal(mu, gene_sigma)
        #Alter gene sizes
        new_gene2_start += round(gene2_eps)
        new_gene2_stop += round(gene2_eps)

        #If gene size is larger or smaller than genome, redo
        while new_gene3_stop > 450 or new_gene1_start < 26:
            gene2_eps = np.random.normal(mu, gene_sigma)
            new_gene2_start = gene2_start_list[-1] + round(gene2_eps)
            new_gene2_stop = gene2_stop_list[-1] + round(gene2_eps)
        while (new_gene2_start > term2_start_list[-1] and new_gene2_start < term2_stop_list[-1]) or (new_gene2_stop > term2_start_list[-1] and new_gene2_stop < term2_stop_list[-1]):
            gene2_eps = np.random.normal(mu, gene_sigma)
            new_gene2_start = gene2_start_list[-1] + round(gene2_eps)
            new_gene2_stop = gene2_stop_list[-1] + round(gene2_eps)
        while (new_gene2_start > prom2_start_list[-1] and new_gene2_start < prom2_stop_list[-1]) or (new_gene2_stop > prom2_start_list[-1] and new_gene2_stop < prom2_stop_list[-1]):
            gene2_eps = np.random.normal(mu, gene_sigma)
            new_gene2_start = gene2_start_list[-1] + round(gene2_eps)
            new_gene2_stop = gene2_stop_list[-1] + round(gene2_eps)
        while (new_gene2_start > rnase2_start_list[-1] and new_gene2_start < rnase2_stop_list[-1]) or (new_gene2_stop > rnase2_start_list[-1] and new_gene2_stop < rnase2_stop_list[-1]):
            gene2_eps = np.random.normal(mu, gene_sigma)
            new_gene2_start = gene2_start_list[-1] + round(gene2_eps)
            new_gene2_stop = gene2_stop_list[-1] + round(gene2_eps)
        #if (new_gene2_start > term2_start_list[-1] and new_gene2_start < term2_stop_list[-1]) or (new_gene2_stop > term2_start_list[-1] and new_gene2_stop < term2_stop_list[-1]) or (new_gene2_start > prom2_start_list[-1] and new_gene2_start < prom2_stop_list[-1]) or (new_gene2_stop > prom2_start_list[-1] and new_gene2_stop < prom2_stop_list[-1]) or (new_gene2_start > rnase2_start_list[-1] and new_gene2_start < rnase2_stop_list[-1]) or (new_gene2_stop > rnase2_start_list[-1] and new_gene2_stop < rnase2_stop_list[-1]):
            #alter_gene2_length(new_gene2_start, new_gene2_stop)

        gene2_start_list.append(new_gene2_start)
        gene2_stop_list.append(new_gene2_stop)

    if chosen_gene == "gene3":
        #Changing gene lengths
        gene3_eps = np.random.normal(mu, gene_sigma)
        #Alter gene sizes
        new_gene3_start += round(gene3_eps)
        new_gene3_stop += round(gene3_eps)
        #If gene size is larger or smaller than genome, redo
        while new_gene3_stop > 450 or new_gene1_start < 26:
            gene3_eps = np.random.normal(mu, gene_sigma)
            new_gene3_start = gene3_start_list[-1] + round(gene3_eps)
            new_gene3_stop = gene3_stop_list[-1] + round(gene3_eps)
        while (new_gene3_start > term3_start_list[-1] and new_gene3_start < term3_stop_list[-1]) or (new_gene3_stop > term3_start_list[-1] and new_gene3_stop < term3_stop_list[-1]):
            gene3_eps = np.random.normal(mu, gene_sigma)
            new_gene3_start = gene3_start_list[-1] + round(gene3_eps)
            new_gene3_stop = gene3_stop_list[-1] + round(gene3_eps)
        while (new_gene3_start > prom2_start_list[-1] and new_gene3_start < prom2_stop_list[-1]):
            gene3_eps = np.random.normal(mu, gene_sigma)
            new_gene3_start = gene3_start_list[-1] + round(gene3_eps)
            new_gene3_stop = gene3_stop_list[-1] + round(gene3_eps)
        while (new_gene3_start > rnase3_start_list[-1] and new_gene3_start < rnase3_stop_list[-1]) or (new_gene3_stop > rnase3_start_list[-1] and new_gene3_stop < rnase3_stop_list[-1]):
            gene3_eps = np.random.normal(mu, gene_sigma)
            new_gene3_start = gene3_start_list[-1] + round(gene3_eps)
            new_gene3_stop = gene3_stop_list[-1] + round(gene3_eps)
        #if (new_gene3_start > term3_start_list[-1] and new_gene3_start < term3_stop_list[-1]) or (new_gene3_stop > term3_start_list[-1] and new_gene3_stop < term3_stop_list[-1]) or (new_gene3_start > prom2_start_list[-1] and new_gene3_start < prom2_stop_list[-1]) or (new_gene3_start > rnase3_start_list[-1] and new_gene3_start < rnase3_stop_list[-1]) or (new_gene3_stop > rnase3_start_list[-1] and new_gene3_stop < rnase3_stop_list[-1]):
            #alter_gene_length3(new_gene3_start, new_gene3_stop)

        gene3_start_list.append(new_gene3_start)
        gene3_stop_list.append(new_gene3_stop)


def add_promoter(new_prom1_start, new_prom1_stop, new_prom2_start, new_prom2_stop):

    prom_sigma = 10.0
    promoter1_start = 0
    promoter2_start = 0
    promoter_possibilities = ["promoter1", "promoter2"]
    chosen_promoter = random.choice(promoter_possibilities)

    if chosen_promoter == "promoter1":
        #Adding in a promoter between genes 1 and 2
        if new_prom1_start > 0:
            promoter1_start = new_prom1_start
            promoter1_stop = new_prom1_stop
        if (new_gene2_start - (new_gene1_stop - 7) >= 25) and (promoter1_start == 0):
            new_prom1_start = random.randint(new_gene1_stop-7, new_gene2_start-25)
            new_prom1_stop = new_prom1_start + 9
            prom1_start_list.append(new_prom1_start)
            prom1_stop_list.append(new_prom1_stop)
        prom1_eps = np.random.normal(mu, prom_sigma)
        #Adjusting promoter site
        new_prom1_start += round(prom1_eps)
        new_prom1_stop += round(prom1_eps)
        prom1_start_list.append(new_prom1_start)
        prom1_stop_list.append(new_prom1_stop)
        while new_prom1_start < new_gene2_start - 15 or new_prom1_start < new_gene1_stop - 7:
            prom1_eps = np.random.normal(mu, prom_sigma)
            new_prom1_start += round(prom1_eps)
            new_prom1_stop += round(prom1_eps)
        while new_prom1_start >= prom2_start_list[-1] and new_prom1_start <= prom2_stop_list[-1]:
            prom1_eps = np.random.normal(mu, prom_sigma)
            new_prom1_start += round(prom1_eps)
            new_prom1_stop += round(prom1_eps)
        while (new_prom1_start > rnase1_start_list[-1] and new_prom1_start < rnase1_stop_list[-1]) or (new_prom1_stop > rnase1_start_list[-1] and new_prom1_stop < rnase1_stop_list[-1]):
            prom1_eps = np.random.normal(mu, prom_sigma)
            new_prom1_start += round(prom1_eps)
            new_prom1_stop += round(prom1_eps)
        while (new_prom1_start > rnase2_start_list[-1] and new_prom1_start < rnase2_stop_list[-1]) or (new_prom1_stop > rnase2_start_list[-1] and new_prom1_stop < rnase2_stop_list[-1]):
            prom1_eps = np.random.normal(mu, prom_sigma)
            new_prom1_start += round(prom1_eps)
            new_prom1_stop += round(prom1_eps)
        while (new_prom1_start > rnase3_start_list[-1] and new_prom1_start < rnase3_stop_list[-1]) or (new_prom1_stop > rnase3_start_list[-1] and new_prom1_stop < rnase3_stop_list[-1]):
            prom1_eps = np.random.normal(mu, prom_sigma)
            new_prom1_start += round(prom1_eps)
            new_prom1_stop += round(prom1_eps)
        while (new_prom1_start > term1_start_list[-1] and new_prom1_start < term1_stop_list[-1]) or (new_prom1_stop > term1_start_list[-1] and new_prom1_stop < term1_stop_list[-1]):
            prom1_eps = np.random.normal(mu, prom_sigma)
            new_prom1_start += round(prom1_eps)
            new_prom1_stop += round(prom1_eps)
        while (new_prom1_start > term2_start_list[-1] and new_prom1_start < term2_stop_list[-1]) or (new_prom1_stop > term2_start_list[-1] and new_prom1_stop < term2_stop_list[-1]):
            prom1_eps = np.random.normal(mu, prom_sigma)
            new_prom1_start += round(prom1_eps)
            new_prom1_stop += round(prom1_eps)
        while (new_prom1_start > term3_start_list[-1] and new_prom1_start < term3_stop_list[-1]) or (new_prom1_stop > term3_start_list[-1] and new_prom1_stop < term3_stop_list[-1]):
            prom1_eps = np.random.normal(mu, prom_sigma)
            new_prom1_start += round(prom1_eps)
            new_prom1_stop += round(prom1_eps)
        promoter1_start = new_prom1_start
        promoter1_stop = new_prom1_stop

    if chosen_promoter == "promoter2":
        #Adding promoter between genes 2 and 3
        if new_prom2_start > 0:
            promoter2_start = new_prom2_start
            promoter2_stop = new_prom2_stop
        if (new_gene3_start - (new_gene2_stop - 7) >= 25) and (promoter2_start == 0):
            new_prom2_start = random.randint(new_gene2_stop-7, new_gene3_start-25)
            new_prom2_stop = new_prom2_start + 9
            prom2_start_list.append(new_prom2_start)
            prom2_stop_list.append(new_prom2_stop)
        prom2_eps = np.random.normal(mu, prom_sigma)
        #Adjusting promoter site
        new_prom2_start += round(prom2_eps)
        new_prom2_stop += round(prom2_eps)
        prom2_start_list.append(new_prom2_start)
        prom2_stop_list.append(new_prom2_stop)
        while new_prom2_start < new_gene3_start - 15 or new_prom2_start < new_gene2_stop - 7:
            prom2_eps = np.random.normal(mu, prom_sigma)
            new_prom2_start += round(prom2_eps)
            new_prom2_stop += round(prom2_eps)
        while new_prom2_start >= prom1_start_list[-1] and new_prom2_start <= prom1_stop_list[-1]:
            prom2_eps = np.random.normal(mu, prom_sigma)
            new_prom2_start += round(prom2_eps)
            new_prom2_stop += round(prom2_eps)
        while (new_prom2_start > rnase1_start_list[-1] and new_prom2_start < rnase1_stop_list[-1]) or (new_prom2_stop > rnase1_start_list[-1] and new_prom2_stop < rnase1_stop_list[-1]):
            prom2_eps = np.random.normal(mu, prom_sigma)
            new_prom2_start += round(prom2_eps)
            new_prom2_stop += round(prom2_eps)
        while (new_prom2_start > rnase2_start_list[-1] and new_prom2_start < rnase2_stop_list[-1]) or (new_prom2_stop > rnase2_start_list[-1] and new_prom2_stop < rnase2_stop_list[-1]):
            prom2_eps = np.random.normal(mu, prom_sigma)
            new_prom2_start += round(prom2_eps)
            new_prom2_stop += round(prom2_eps)
        while (new_prom2_start > rnase3_start_list[-1] and new_prom2_start < rnase3_stop_list[-1]) or (new_prom2_stop > rnase2_start_list[-1] and new_prom2_stop < rnase3_stop_list[-1]):
            prom2_eps = np.random.normal(mu, prom_sigma)
            new_prom2_start += round(prom2_eps)
            new_prom2_stop += round(prom2_eps)
        while (new_prom2_start > term1_start_list[-1] and new_prom2_start < term1_stop_list[-1]) or (new_prom2_stop > term1_start_list[-1] and new_prom2_stop < term1_stop_list[-1]):
            prom2_eps = np.random.normal(mu, prom_sigma)
            new_prom2_start += round(prom1_eps)
            new_prom2_stop += round(prom1_eps)
        while (new_prom2_start > term2_start_list[-1] and new_prom2_start < term2_stop_list[-1]) or (new_prom2_stop > term2_start_list[-1] and new_prom2_stop < term2_stop_list[-1]):
            prom2_eps = np.random.normal(mu, prom_sigma)
            new_prom2_start += round(prom1_eps)
            new_prom2_stop += round(prom1_eps)
        while (new_prom2_start > term3_start_list[-1] and new_prom2_start < term3_stop_list[-1]) or (new_prom2_stop > term3_start_list[-1] and new_prom2_stop < term3_stop_list[-1]):
            prom2_eps = np.random.normal(mu, prom_sigma)
            new_prom2_start += round(prom1_eps)
            new_prom2_stop += round(prom1_eps)
        promoter2_start = new_prom2_start
        promoter2_stop = new_prom2_stop

def remove_promoter(promoter1_start, promoter1_stop, promoter2_start, promoter2_stop):

    promoter_possibilities = ["promoter1", "promoter2"]
    chosen_promoter = random.choice(promoter_possibilities)

    if chosen_promoter == "promoter1":
        #Removing promoter between genes 1 and 2
        promoter1_start = 0
        promoter1_stop = 0

    if chosen_promoter == "promoter2":
        #Removing promoter between genes 2 and 3
        promoter2_start = 0
        promoter2_stop = 0

def add_ranse(promoter1_start, promoter2_start):

    rnase_possibilities = ["rnase1", "rnase2", "rnase3"]
    chosen_rnase = random.choice(rnase_possibilities)

    if chosen_rnase == "rnase1":
        #Adds rnase after first promoter
        new_rnase1_start = 10
        new_rnase1_stop = 20
        if (new_rnase1_start > term1_start_list[-1] and new_rnase1_start < term1_stop_list[-1]) or (new_rnase1_stop > term1_start_list[-1] and new_rnase1_stop < term1_stop_list[-1]):
            new_rnase1_start = 0
            new_rnase1_stop = 0
        if (new_rnase1_start > term2_start_list[-1] and new_rnase1_start < term2_stop_list[-1]) or (new_rnase1_stop > term2_start_list[-1] and new_rnase1_stop < term2_stop_list[-1]):
            new_rnase1_start = 0
            new_rnase1_stop = 0
        if (new_rnase1_start > term3_start_list[-1] and new_rnase1_start < term3_stop_list[-1]) or (new_rnase1_stop > term3_start_list[-1] and new_rnase1_stop < term3_stop_list[-1]):
            new_rnase1_start = 0
            new_rnase1_stop = 0
        rnase1_start_list.append(new_rnase1_start)
        rnase1_stop_list.append(new_rnase1_stop)

    if chosen_rnase == "rnase2":
        #Adds rnase after second promoter
        new_rnase2_start = 0
        new_rnase2_stop = 0
        if promoter1_start > 0:
            new_rnase2_start = promoter1_stop
            new_rnase2_stop = new_rnase2_start + 10
        if (new_rnase2_start > prom1_start_list[-1] and new_rnase2_start < prom1_stop_list[-1]) or (new_rnase2_stop > prom1_start_list[-1] and new_rnase2_stop < prom1_stop_list[-1]):
            new_rnase2_start = 0
            new_rnase2_stop = 0
        if (new_rnase2_start > prom2_start_list[-1] and new_rnase2_start < prom2_stop_list[-1]) or (new_rnase2_stop > prom2_start_list[-1] and new_rnase2_stop < prom2_stop_list[-1]):
            new_rnase2_start = 0
            new_rnase2_stop = 0
        if (new_rnase2_start > term1_start_list[-1] and new_rnase2_start < term1_stop_list[-1]) or (new_rnase2_stop > term1_start_list[-1] and new_rnase2_stop < term1_stop_list[-1]):
            new_rnase2_start = 0
            new_rnase2_stop = 0
        if (new_rnase2_start > term2_start_list[-1] and new_rnase2_start < term2_stop_list[-1]) or (new_rnase2_stop > term2_start_list[-1] and new_rnase2_stop < term2_stop_list[-1]):
            new_rnase2_start = 0
            new_rnase2_stop = 0
        if (new_rnase2_start > term3_start_list[-1] and new_rnase2_start < term3_stop_list[-1]) or (new_rnase2_stop > term3_start_list[-1] and new_rnase2_stop < term3_stop_list[-1]):
            new_rnase2_start = 0
            new_rnase2_stop = 0
        rnase2_start_list.append(new_rnase2_start)
        rnase2_stop_list.append(new_rnase2_stop)

    if chosen_rnase == "rnase3":
        #Adds rnase after third promoter
        new_rnase3_start = 0
        new_rnase3_stop = 0
        if promoter2_start > 0:
            new_rnase3_start = promoter2_stop
            new_rnase3_stop = new_rnase3_start + 10
        if (new_rnase3_start > prom1_start_list[-1] and new_rnase3_start < prom1_stop_list[-1]) or (new_rnase3_stop > prom1_start_list[-1] and new_rnase3_stop < prom1_stop_list[-1]):
            new_rnase3_start = 0
            new_rnase3_stop = 0
        if (new_rnase3_start > prom2_start_list[-1] and new_rnase3_start < prom2_stop_list[-1]) or (new_rnase3_stop > prom2_start_list[-1] and new_rnase3_stop < prom2_stop_list[-1]):
            new_rnase3_start = 0
            new_rnase3_stop = 0
        if (new_rnase3_start > term1_start_list[-1] and new_rnase3_start < term1_stop_list[-1]) or (new_rnase3_stop > term1_start_list[-1] and new_rnase3_stop < term1_stop_list[-1]):
            new_rnase3_start = 0
            new_rnase3_stop = 0
        if (new_rnase3_start > term2_start_list[-1] and new_rnase3_start < term2_stop_list[-1]) or (new_rnase3_stop > term2_start_list[-1] and new_rnase3_stop < term2_stop_list[-1]):
            new_rnase3_start = 0
            new_rnase3_stop = 0
        if (new_rnase3_start > term3_start_list[-1] and new_rnase3_start < term3_stop_list[-1]) or (new_rnase3_stop > term3_start_list[-1] and new_rnase3_stop < term3_stop_list[-1]):
            new_rnase3_start = 0
            new_rnase3_stop = 0
        rnase3_start_list.append(new_rnase3_start)
        rnase3_stop_list.append(new_rnase3_stop)

def add_terminator(new_gene2_start, new_gene1_stop, new_gene3_start, new_gene2_stop, new_gene3_stop):

    terminator_possibilities = ["terminator1", "terminator2", "terminator3"]
    chosen_terminator = random.choice(terminator_possibilities)

    if chosen_terminator == "terminator1":
        #Adds terminator after first gene
        new_term1_start = 0
        new_term1_stop = 0
        if new_gene2_start - new_gene1_stop > 16:
            new_term1_start = new_gene1_stop + 1
            new_term1_stop = new_gene1_stop + 2
        if (new_term1_start > prom1_start_list[-1] and new_term1_start < prom1_stop_list[-1]) or (new_term1_stop > prom1_start_list[-1] and new_term1_stop < prom1_stop_list[-1]):
            new_term1_start = 0
            new_term1_stop = 0
        if (new_term1_start > prom2_start_list[-1] and new_term1_start < prom2_stop_list[-1]) or (new_term1_stop > prom2_start_list[-1] and new_term1_stop < prom2_stop_list[-1]):
            new_term1_start = 0
            new_term1_stop = 0
        if (new_term1_start > rnase1_start_list[-1] and new_term1_start < rnase1_stop_list[-1]) or (new_term1_stop > rnase1_start_list[-1] and new_term1_stop < rnase1_stop_list[-1]):
            new_term1_start = 0
            new_term1_stop = 0
        if (new_term1_start > rnase2_start_list[-1] and new_term1_start < rnase2_stop_list[-1]) or (new_term1_stop > rnase2_start_list[-1] and new_term1_stop < rnase2_stop_list[-1]):
            new_term1_start = 0
            new_term1_stop = 0
        if (new_term1_start > rnase3_start_list[-1] and new_term1_start < rnase3_stop_list[-1]) or (new_term1_stop > rnase2_start_list[-1] and new_term1_stop < rnase3_stop_list[-1]):
            new_term1_start = 0
            new_term1_stop = 0
        term1_start_list.append(new_term1_start)
        term1_stop_list.append(new_term1_stop)

    if chosen_terminator == "terminator2":
        #Adds terminator after second gene
        new_term2_start = 0
        new_term2_stop = 0
        if new_gene3_start - new_gene2_stop > 16:
            new_term2_start = new_gene2_stop + 1
            new_term2_stop = new_gene2_stop + 2
        if (new_term2_start > prom1_start_list[-1] and new_term2_start < prom1_stop_list[-1]) or (new_term2_stop > prom1_start_list[-1] and new_term2_stop < prom1_stop_list[-1]):
            new_term2_start = 0
            new_term2_stop = 0
        if (new_term2_start > prom2_start_list[-1] and new_term2_start < prom2_stop_list[-1]) or (new_term2_stop > prom2_start_list[-1] and new_term2_stop < prom2_stop_list[-1]):
            new_term2_start = 0
            new_term2_stop = 0
        if (new_term2_start > rnase1_start_list[-1] and new_term2_start < rnase1_stop_list[-1]) or (new_term2_stop > rnase1_start_list[-1] and new_term2_stop < rnase1_stop_list[-1]):
            new_term2_start = 0
            new_term2_stop = 0
        if (new_term2_start > rnase2_start_list[-1] and new_term2_start < rnase2_stop_list[-1]) or (new_term2_stop > rnase2_start_list[-1] and new_term2_stop < rnase2_stop_list[-1]):
            new_term2_start = 0
            new_term2_stop = 0
        if (new_term2_start > rnase3_start_list[-1] and new_term2_start < rnase3_stop_list[-1]) or (new_term2_stop > rnase2_start_list[-1] and new_term2_stop < rnase3_stop_list[-1]):
            new_term2_start = 0
            new_term2_stop = 0
        term2_start_list.append(new_term2_start)
        term2_stop_list.append(new_term2_stop)

    if chosen_terminator == "terminator3":
        #Adds terminator after third gene
        new_term3_start = 0
        new_term3_stop = 0
        if new_gene3_stop < 449:
            new_term3_start = new_gene3_stop + 1
            new_term3_stop = new_gene3_stop + 2
        if (new_term3_start > prom1_start_list[-1] and new_term3_start < prom1_stop_list[-1]) or (new_term3_stop > prom1_start_list[-1] and new_term3_stop < prom1_stop_list[-1]):
            new_term3_start = 0
            new_term3_stop = 0
        if (new_term3_start > prom2_start_list[-1] and new_term3_start < prom2_stop_list[-1]) or (new_term3_stop > prom2_start_list[-1] and new_term3_stop < prom2_stop_list[-1]):
            new_term3_start = 0
            new_term3_stop = 0
        if (new_term3_start > rnase1_start_list[-1] and new_term3_start < rnase1_stop_list[-1]) or (new_term3_stop > rnase1_start_list[-1] and new_term3_stop < rnase1_stop_list[-1]):
            new_term3_start = 0
            new_term3_stop = 0
        if (new_term3_start > rnase2_start_list[-1] and new_term3_start < rnase2_stop_list[-1]) or (new_term3_stop > rnase2_start_list[-1] and new_term3_stop < rnase2_stop_list[-1]):
            new_term3_start = 0
            new_term3_stop = 0
        if (new_term3_start > rnase3_start_list[-1] and new_term3_start < rnase3_stop_list[-1]) or (new_term3_stop > rnase2_start_list[-1] and new_term3_stop < rnase2_stop_list[-1]):
            new_term3_start = 0
            new_term3_stop = 0
        term3_start_list.append(new_term3_start)
        term3_stop_list.append(new_term3_stop)

def remove_terminator(new_term1_start, new_term1_stop, new_term2_start, new_term2_stop, new_term3_start, new_term3_stop):

    terminator_possibilities = ["terminator1", "terminator2", "terminator3"]
    chosen_terminator = random.choice(terminator_possibilities)

    if chosen_terminator == "terminator1":
        #Removes terminator after first gene
        new_term1_start = 0
        new_term1_stop = 0

    if chosen_terminator == "terminator2":
        #Removes terminator after second gene
        new_term2_start = 0
        new_term2_stop = 0

    if chosen_terminator == "terminator3":
        #Removes terminator after third gene
        new_term3_start = 0
        new_term3_stop = 0

def remove_rnase(new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop, new_rnase3_start, new_rnase3_stop):

    rnase_possibilities = ["rnase1", "rnase2", "rnase3"]
    chosen_rnase = random.choice(rnase_possibilities)

    if chosen_rnase == "rnase1":
        #Removes rnase after first promoter
        new_rnase1_start = 0
        new_rnase1_stop = 0

    if chosen_rnase == "rnase2":
        #Removes rnase after second promoter
        new_rnase2_start = 0
        new_rnase2_stop = 0

    if chosen_rnase == "rnase3":
        #Removes rnase after third promoter
        new_rnase3_start = 0
        new_rnase3_stop = 0

#Class containing genome for simulation with new polymerase strength
class three_genome:

    def __init__(self):
        self.pol_strength = 4e10

    def recreated_genome(pol_strength, term_efficiency, gene1_start, gene1_stop, gene2_start,
                         gene2_stop, gene3_start, gene3_stop, prom1_start, prom1_stop, prom2_start,
                         prom2_stop, rnase1_start, rnase1_stop, rnase2_start, rnase2_stop, rnase3_start,
                         rnase3_stop, term1_start, term1_stop, term2_start, term2_stop, term3_start,
                         term3_stop):

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
                                 interactions={"rnapol": pol_strength})
        if prom2_start > 0:
            plasmid.add_promoter(name="p3", start=prom2_start, stop=prom2_stop,
                                 interactions={"rnapol": pol_strength})
        if rnase1_start > 0:
            plasmid.add_rnase_site(start=rnase1_start, stop=rnase1_stop)
        if rnase2_start > 0:
            plasmid.add_rnase_site(start=rnase2_start, stop=rnase2_stop)
        if rnase3_start > 0:
            plasmid.add_rnase_site(start=rnase3_start, stop=rnase3_stop)
        if term1_start > 0:
            plasmid.add_terminator(name="t1", start=term1_start, stop=term1_stop,
                                   efficiency={"rnapol": term_efficiency})
        if term2_start > 0:
            plasmid.add_terminator(name="t2", start=term2_start, stop=term2_stop,
                                   efficiency={"rnapol": term_efficiency})
        if term3_start > 0:
            plasmid.add_terminator(name="t3", start=term3_start, stop=term3_stop,
                                   efficiency={"rnapol": term_efficiency})
        plasmid.add_gene(name="proteinX", start=gene1_start, stop=gene1_stop,
                         rbs_start=gene1_start-15, rbs_stop=gene1_start, rbs_strength=1e7)

        plasmid.add_gene(name="proteinY", start=gene2_start, stop=gene2_stop,
                         rbs_start=gene2_start-15, rbs_stop=gene2_start, rbs_strength=1e7)
        plasmid.add_gene(name="proteinZ", start=gene3_start, stop=gene3_stop,
                         rbs_start=gene3_start-15, rbs_stop=gene3_start, rbs_strength=1e7)
        sim.register_genome(plasmid)
        sim.simulate(time_limit=240, time_step=1,
                     output = "three_genes_replicated.tsv")

    def best_recreated_genome(pol_strength, term_efficiency, gene1_start, gene1_stop, gene2_start,
                              gene2_stop, gene3_start, gene3_stop, prom1_start, prom1_stop, prom2_start,
                              prom2_stop, rnase1_start, rnase1_stop, rnase2_start, rnase2_stop, rnase3_start,
                              rnase3_stop, term1_start, term1_stop, term2_start, term2_stop, term3_start,
                              term3_stop):

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
                                 interactions={"rnapol": pol_strength})
        if prom2_start > 0:
            plasmid.add_promoter(name="p3", start=prom2_start, stop=prom2_stop,
                                 interactions={"rnapol": pol_strength})
        if rnase1_start > 0:
            plasmid.add_rnase_site(start=rnase1_start, stop=rnase1_stop)
        if rnase2_start > 0:
            plasmid.add_rnase_site(start=rnase2_start, stop=rnase2_stop)
        if rnase3_start > 0:
            plasmid.add_rnase_site(start=rnase3_start, stop=rnase3_stop)
        if term1_start > 0:
            plasmid.add_terminator(name="t1", start=term1_start, stop=term1_stop,
                                   efficiency={"rnapol": term_efficiency})
        if term2_start > 0:
            plasmid.add_terminator(name="t2", start=term2_start, stop=term2_stop,
                                   efficiency={"rnapol": term_efficiency})
        if term3_start > 0:
            plasmid.add_terminator(name="t3", start=term3_start, stop=term3_stop,
                                   efficiency={"rnapol": term_efficiency})
        plasmid.add_gene(name="proteinX", start=gene1_start, stop=gene1_stop,
                         rbs_start=gene1_start-15, rbs_stop=gene1_start, rbs_strength=1e7)

        plasmid.add_gene(name="proteinY", start=gene2_start, stop=gene2_stop,
                         rbs_start=gene2_start-15, rbs_stop=gene2_start, rbs_strength=1e7)
        plasmid.add_gene(name="proteinZ", start=gene3_start, stop=gene3_stop,
                         rbs_start=gene3_start-15, rbs_stop=gene3_start, rbs_strength=1e7)
        sim.register_genome(plasmid)
        sim.simulate(time_limit=240, time_step=1,
                     output = "best_three_genes_replicated.tsv")


if __name__ == '__main__':
    main()
