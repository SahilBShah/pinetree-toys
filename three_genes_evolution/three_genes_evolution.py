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
    prom1_start_list = [0]
    prom1_stop_list = [0]
    prom2_start_list = [0]
    prom2_stop_list = [0]
    f_old = 1.0
    mu = 0.0
    sigma = 1.0
    Ne = 10
    i = 0
    gen = 0
    new_gene1_start = 26
    new_gene1_stop = 121
    '''new_gene2_start = 178
    new_gene2_stop = 300
    new_gene3_start = 315
    new_gene3_stop = 448'''
    new_gene2_start = 159
    new_gene2_stop = 280
    new_gene3_start = 319
    new_gene3_stop = 449
    '''promoter1_start = 0
    promoter1_stop = 0
    promoter2_start = 0
    promoter2_stop = 0'''
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
    '''new_gene1_start = random.randint(25, 449)
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
    gene3_start_list.append(new_gene3_start)'''

    '''possibilities = ["alter polymerase strength", "alter terminator efficiency", "alter gene length", "add promoter", "remove promoter", "add rnase",
                        "remove rnase", "add terminator", "remove terminator"]'''

    possibilities = ["alter polymerase strength", "alter terminator efficiency", "add promoter", "remove promoter", "add rnase",
                      "remove rnase", "add terminator", "remove terminator"]

    while i < 10000:

        mutation = random.choice(possibilities)

        if mutation == "alter polymerase strength":
            new_pol_strength = alter_poly_strength(new_pol_strength)
        if mutation == "alter terminator efficiency":
            new_term_efficiency = alter_term_efficiency(new_term_efficiency)
        if mutation == "alter gene length":
            gene_length = alter_gene_length(new_gene1_start, new_gene1_stop, new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop)
            new_gene1_start = gene1_start_list[-1]
            new_gene1_stop = gene1_stop_list[-1]
            new_gene2_start = gene2_start_list[-1]
            new_gene2_stop = gene2_stop_list[-1]
            new_gene3_start = gene3_start_list[-1]
            new_gene3_stop = gene3_stop_list[-1]
        if mutation == "add promoter":
            new_promoter = add_promoter(new_prom1_start, new_prom1_stop, new_prom2_start, new_prom2_stop)
            new_prom1_start = prom1_start_list[-1]
            new_prom1_stop = prom1_stop_list[-1]
            new_prom2_start = prom2_start_list[-1]
            new_prom2_stop = prom2_stop_list[-1]
        if mutation == "remove promoter":
            old_promoter = remove_promoter(new_prom1_start, new_prom1_stop, new_prom2_start, new_prom2_stop)
            new_prom1_start = prom1_start_list[-1]
            new_prom1_stop = prom1_stop_list[-1]
            new_prom2_start = prom2_start_list[-1]
            new_prom2_stop = prom2_stop_list[-1]
        if mutation == "add rnase":
            new_rnase = add_ranse(new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop, new_rnase3_start, new_rnase3_stop)
            new_rnase1_start = rnase1_start_list[-1]
            new_rnase1_stop = rnase1_stop_list[-1]
            new_rnase2_start = rnase2_start_list[-1]
            new_rnase2_stop = rnase2_stop_list[-1]
            new_rnase3_start = rnase3_start_list[-1]
            new_rnase3_stop = rnase3_stop_list[-1]
        if mutation == "remove rnase":
            old_rnase = remove_rnase(new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop, new_rnase3_start, new_rnase3_stop)
            new_rnase1_start = rnase1_start_list[-1]
            new_rnase1_stop = rnase1_stop_list[-1]
            new_rnase2_start = rnase2_start_list[-1]
            new_rnase2_stop = rnase2_stop_list[-1]
            new_rnase3_start = rnase3_start_list[-1]
            new_rnase3_stop = rnase3_stop_list[-1]
        if mutation == "add terminator":
            new_terminator = add_terminator(new_term1_start, new_term1_stop, new_term2_start, new_term2_stop, new_term3_start, new_term3_stop)
            new_term1_start = term1_start_list[-1]
            new_term1_stop = term1_stop_list[-1]
            new_term2_start = term2_start_list[-1]
            new_term2_stop = term2_stop_list[-1]
            new_term3_start = term3_start_list[-1]
            new_term3_stop = term3_stop_list[-1]
        if mutation == "remove terminator":
            old_terminator = remove_terminator(new_term1_start, new_term1_stop, new_term2_start, new_term2_stop, new_term3_start, new_term3_stop)
            new_term1_start = term1_start_list[-1]
            new_term1_stop = term1_stop_list[-1]
            new_term2_start = term2_start_list[-1]
            new_term2_stop = term2_stop_list[-1]
            new_term3_start = term3_start_list[-1]
            new_term3_stop = term3_stop_list[-1]

        eps = np.random.normal(mu, sigma)
        f_new = f_old * (1.0 + eps)
        test = accept_mutation(df, f_old, f_new, new_pol_strength, new_term_efficiency, new_gene1_start, new_gene1_stop, new_gene2_start,
                                    new_gene2_stop, new_gene3_start, new_gene3_stop, new_prom1_start, new_prom1_stop, new_prom2_start,
                                    new_prom2_stop, new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop, new_rnase3_start,
                                    new_rnase3_stop, new_term1_start, new_term1_stop, new_term2_start, new_term2_stop, new_term3_start,
                                    new_term3_stop)
        if test > 0:
            f_old = test
        '''if sos_list[-1] == 0:
            break'''

        i+=1
        print("i = ", i)


    #Exported tsv files
    all_sos_dataframe = pandas.DataFrame(all_sos_list, columns=["Sum_of_Squares"])
    export_csv = all_sos_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/all_sos_data.tsv", index=False)
    sos_dataframe = pandas.DataFrame(sos_list, columns=["Sum_of_Squares"])
    export_csv = sos_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/sos_data.tsv", index=False)
    all_poly_dataframe = pandas.DataFrame(all_poly_list, columns=["Polymerase_Rate"])
    export_csv = all_poly_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/all_poly_data.tsv", index=False)
    poly_dataframe = pandas.DataFrame(poly_list, columns=["Polymerase_Rate"])
    export_csv = poly_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/poly_data.tsv", index=False)
    all_term_dataframe = pandas.DataFrame(all_term_list, columns=["Terminator Efficiency Rate"])
    export_csv = all_term_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/all_term_data.tsv", index=False)
    term_dataframe = pandas.DataFrame(term_list, columns=["Terminator_Efficiency_Rate"])
    export_csv = term_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/term_data.tsv", index=False)
    gene1_start_dataframe = pandas.DataFrame(gene1_start_list, columns=["Gene_Start"])
    export_csv = gene1_start_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/gene1_start_data.tsv", index=False)
    gene2_start_dataframe = pandas.DataFrame(gene2_start_list, columns=["Gene_Start"])
    export_csv = gene2_start_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/gene2_start_data.tsv", index=False)
    gene3_start_dataframe = pandas.DataFrame(gene3_start_list, columns=["Gene Start"])
    export_csv = gene3_start_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/gene3_start_data.tsv", index=False)
    gene1_stop_dataframe = pandas.DataFrame(gene1_stop_list, columns=["Gene_Stop"])
    export_csv = gene1_stop_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/gene1_stop_data.tsv", index=False)
    gene2_stop_dataframe = pandas.DataFrame(gene2_stop_list, columns=["Gene_Stop"])
    export_csv = gene2_stop_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/gene2_stop_data.tsv", index=False)
    gene3_stop_dataframe = pandas.DataFrame(gene3_stop_list, columns=["Gene_Stop"])
    export_csv = gene3_stop_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/gene3_stop_data.tsv", index=False)
    print("Final Polymerase Rate: ", poly_list[-1])
    print("Final Terminator Efficiency Rate: ", term_list[-1])
    print("Gene 1 starts at: ", new_gene1_start, "and stops at: ", new_gene1_stop)
    print("Gene 1 Ribosome Binding Site starts at: ", new_gene1_start-15, "and stops at: ", new_gene1_start)
    print("Gene 1 terminator starts at: ", new_term1_start, "and stops at: ", new_term1_stop)
    print("Gene 2 starts at: ", new_gene2_start, "and stops at: ", new_gene2_stop)
    print("Gene 2 Ribosome Binding Site starts at: ", new_gene2_start-15, "and stops at: ", new_gene2_start)
    print("Gene 2 terminator starts at: ", new_term2_start, "and stops at: ", new_term2_stop)
    print("Gene 3 starts at: ", new_gene3_start, "and stops at: ", new_gene3_stop)
    print("Gene 3 Ribosome Binding Site starts at: ", new_gene3_start-15, "and stops at: ", new_gene3_start)
    print("Gene 3 terminator starts at: ", term3_start_list[-1], "and stops at: ", term3_stop_list[-1])
    print("Promoter 1 starts at: ", 1, "and stops at: ", 10)
    print("Added Promoter 2 starts at: ", prom1_start_list[-1], "and stops at: ", prom1_stop_list[-1])
    print("Added Promoter 3 starts at: ", prom2_start_list[-1], "and stops at: ", prom2_stop_list[-1])
    print("Rnase site 1 starts at: ", rnase3_start_list[-1], "and stops at: ", rnase3_stop_list[-1])
    print("Rnase site 2 starts at: ", rnase2_start_list[-1], "and stops at: ", rnase2_stop_list[-1])
    print("Rnase site 3 starts at: ", rnase1_start_list[-1], "and stops at: ", rnase1_stop_list[-1])
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
        p = ((1-pow((xi/xj), 2)) / (1-pow((xi/xj), (2 * float(N)))))
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

def alter_poly_strength(pol_strength):

    poly_sigma = 1e10
    #Determining polymerase strength
    poly_eps = np.random.normal(mu, poly_sigma)
    #Determines new polymerase strength to reduce sum of squares value
    pol_strength = poly_list[-1] + poly_eps
    while pol_strength < 0:
        poly_eps = np.random.normal(mu, poly_sigma)
        pol_strength = poly_list[-1] + poly_eps
        if pol_strength in all_poly_list:
            poly_eps = np.random.normal(mu, poly_sigma)
            pol_strength = poly_list[-1] + poly_eps
    if pol_strength in all_poly_list:
        poly_eps = np.random.normal(mu, poly_sigma)
        pol_strength = poly_list[-1] + poly_eps
    all_poly_list.append(pol_strength)
    return pol_strength

def alter_term_efficiency(term_efficiency):

    term_sigma = 1.0
    #Determining terminator polymerase efficiency rate
    term_eps = np.random.normal(mu, term_sigma)
    #Determines new terminator efficiency value to reduce sum of squares value
    term_efficiency = term_list[-1] + term_eps
    while term_efficiency < 0 or term_efficiency > 1:
        term_eps = np.random.normal(mu, term_sigma)
        term_efficiency = term_list[-1] + term_eps
        if term_efficiency in all_term_list:
            term_eps = np.random.normal(mu, term_sigma)
            term_efficiency = term_list[-1] + term_eps
    if term_efficiency in all_term_list:
        term_eps = np.random.normal(mu, term_sigma)
        term_efficiency = term_list[-1] + term_eps
    all_term_list.append(term_efficiency)
    return term_efficiency

def alter_gene_length(gene1_start, gene1_stop, gene2_start, gene2_stop, gene3_start, gene3_stop):

    gene_sigma = 10.0
    gene_possibilities = ["gene1", "gene2", "gene3"]
    chosen_gene = random.choice(gene_possibilities)
    if chosen_gene == "gene1":
        #Changing gene lengths
        gene1_eps = np.random.normal(mu, gene_sigma)
        #Alter gene sizes
        gene1_start += round(gene1_eps)
        gene1_stop += round(gene1_eps)

        #If gene size is larger or smaller than genome, redo
        while gene3_stop > 450 or gene1_start < 26:
            #gene1_eps = np.random.normal(mu, gene_sigma)
            gene1_start = gene1_start_list[-1] + round(gene1_eps)
            gene1_stop = gene1_stop_list[-1] + round(gene1_eps)
        while (gene1_start > term1_start_list[-1] and gene1_start < term1_stop_list[-1]) or (gene1_stop > term1_start_list[-1] and gene1_stop < term1_stop_list[-1]):
            #gene1_eps = np.random.normal(mu, gene_sigma)
            gene1_start = gene1_start_list[-1] + round(gene1_eps)
            gene1_stop = gene1_stop_list[-1] + round(gene1_eps)
        while (gene1_start > prom1_start_list[-1] and gene1_start < prom1_stop_list[-1]) or (gene1_stop > prom1_start_list[-1] and gene1_stop < prom1_stop_list[-1]):
            #gene1_eps = np.random.normal(mu, gene_sigma)
            gene1_start = gene1_start_list[-1] + round(gene1_eps)
            gene1_stop = gene1_stop_list[-1] + round(gene1_eps)
        while (gene1_start > rnase1_start_list[-1] and gene1_start < rnase1_stop_list[-1]) or (gene1_stop > rnase1_start_list[-1] and gene1_stop < rnase1_stop_list[-1]):
            #gene1_eps = np.random.normal(mu, gene_sigma)
            gene1_start = gene1_start_list[-1] + round(gene1_eps)
            gene1_stop = gene1_stop_list[-1] + round(gene1_eps)
        #if (new_gene1_start > prom1_start_list[-1] and new_gene1_start < prom1_stop_list[-1]) or (new_gene1_stop > prom1_start_list[-1] and new_gene1_stop < prom1_stop_list[-1]) or (new_gene1_start > prom1_start_list[-1] and new_gene1_start < prom1_stop_list[-1]) or (new_gene1_stop > prom1_start_list[-1] and new_gene1_stop < prom1_stop_list[-1]) or (new_gene1_start > rnase1_start_list[-1] and new_gene1_start < rnase1_stop_list[-1]) or (new_gene1_stop > rnase1_start_list[-1] and new_gene1_stop < rnase1_stop_list[-1]) or new_gene3_stop > 450 or new_gene1_start < 26:
            #alter_gene1_length(new_gene1_start, new_gene1_stop)

        gene1_start_list.append(gene1_start)
        gene1_stop_list.append(gene1_stop)

    if chosen_gene == "gene2":
        #Changing gene lengths
        gene2_eps = np.random.normal(mu, gene_sigma)
        #Alter gene sizes
        gene2_start += round(gene2_eps)
        gene2_stop += round(gene2_eps)

        #If gene size is larger or smaller than genome, redo
        while gene3_stop > 450 or gene1_start < 26:
            #gene2_eps = np.random.normal(mu, gene_sigma)
            gene2_start = gene2_start_list[-1] + round(gene2_eps)
            gene2_stop = gene2_stop_list[-1] + round(gene2_eps)
        while (gene2_start > term2_start_list[-1] and gene2_start < term2_stop_list[-1]) or (gene2_stop > term2_start_list[-1] and gene2_stop < term2_stop_list[-1]):
            #gene2_eps = np.random.normal(mu, gene_sigma)
            gene2_start = gene2_start_list[-1] + round(gene2_eps)
            gene2_stop = gene2_stop_list[-1] + round(gene2_eps)
        while (gene2_start > prom2_start_list[-1] and gene2_start < prom2_stop_list[-1]) or (gene2_stop > prom2_start_list[-1] and gene2_stop < prom2_stop_list[-1]):
            #gene2_eps = np.random.normal(mu, gene_sigma)
            gene2_start = gene2_start_list[-1] + round(gene2_eps)
            gene2_stop = gene2_stop_list[-1] + round(gene2_eps)
        while (gene2_start > rnase2_start_list[-1] and gene2_start < rnase2_stop_list[-1]) or (gene2_stop > rnase2_start_list[-1] and gene2_stop < rnase2_stop_list[-1]):
            #gene2_eps = np.random.normal(mu, gene_sigma)
            gene2_start = gene2_start_list[-1] + round(gene2_eps)
            gene2_stop = gene2_stop_list[-1] + round(gene2_eps)
        #if (new_gene2_start > term2_start_list[-1] and new_gene2_start < term2_stop_list[-1]) or (new_gene2_stop > term2_start_list[-1] and new_gene2_stop < term2_stop_list[-1]) or (new_gene2_start > prom2_start_list[-1] and new_gene2_start < prom2_stop_list[-1]) or (new_gene2_stop > prom2_start_list[-1] and new_gene2_stop < prom2_stop_list[-1]) or (new_gene2_start > rnase2_start_list[-1] and new_gene2_start < rnase2_stop_list[-1]) or (new_gene2_stop > rnase2_start_list[-1] and new_gene2_stop < rnase2_stop_list[-1]):
            #alter_gene2_length(new_gene2_start, new_gene2_stop)

        gene2_start_list.append(gene2_start)
        gene2_stop_list.append(gene2_stop)

    if chosen_gene == "gene3":
        #Changing gene lengths
        gene3_eps = np.random.normal(mu, gene_sigma)
        #Alter gene sizes
        gene3_start += round(gene3_eps)
        gene3_stop += round(gene3_eps)
        #If gene size is larger or smaller than genome, redo
        while gene3_stop > 450 or gene1_start < 26:
            #gene3_eps = np.random.normal(mu, gene_sigma)
            gene3_start = gene3_start_list[-1] + round(gene3_eps)
            gene3_stop = gene3_stop_list[-1] + round(gene3_eps)
        while (gene3_start > term3_start_list[-1] and gene3_start < term3_stop_list[-1]) or (gene3_stop > term3_start_list[-1] and gene3_stop < term3_stop_list[-1]):
            #gene3_eps = np.random.normal(mu, gene_sigma)
            gene3_start = gene3_start_list[-1] + round(gene3_eps)
            gene3_stop = gene3_stop_list[-1] + round(gene3_eps)
        while (gene3_start > prom2_start_list[-1] and gene3_start < prom2_stop_list[-1]):
            #gene3_eps = np.random.normal(mu, gene_sigma)
            gene3_start = gene3_start_list[-1] + round(gene3_eps)
            gene3_stop = gene3_stop_list[-1] + round(gene3_eps)
        while (gene3_start > rnase3_start_list[-1] and gene3_start < rnase3_stop_list[-1]) or (gene3_stop > rnase3_start_list[-1] and gene3_stop < rnase3_stop_list[-1]):
            #gene3_eps = np.random.normal(mu, gene_sigma)
            gene3_start = gene3_start_list[-1] + round(gene3_eps)
            gene3_stop = gene3_stop_list[-1] + round(gene3_eps)
        #if (new_gene3_start > term3_start_list[-1] and new_gene3_start < term3_stop_list[-1]) or (new_gene3_stop > term3_start_list[-1] and new_gene3_stop < term3_stop_list[-1]) or (new_gene3_start > prom2_start_list[-1] and new_gene3_start < prom2_stop_list[-1]) or (new_gene3_start > rnase3_start_list[-1] and new_gene3_start < rnase3_stop_list[-1]) or (new_gene3_stop > rnase3_start_list[-1] and new_gene3_stop < rnase3_stop_list[-1]):
            #alter_gene_length3(new_gene3_start, new_gene3_stop)

        gene3_start_list.append(gene3_start)
        gene3_stop_list.append(gene3_stop)


def add_promoter(prom1_start, prom1_stop, prom2_start, prom2_stop):

    region1a_start = 122
    region1a_stop = 132
    region1b_start = 133
    region1b_stop = 143
    region2a_start = 281
    region2a_stop = 291
    region2b_start = 292
    region2b_stop = 302
    promoter_possibilities = ["promoter1", "promoter2"]
    chosen_promoter = random.choice(promoter_possibilities)
    promoter1_slots = ['A', 'B']
    promoter2_slots = ['A', 'B']

    if chosen_promoter == "promoter1":
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

    if chosen_promoter == "promoter2":
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

def remove_promoter(promoter1_start, promoter1_stop, promoter2_start, promoter2_stop):

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

def add_ranse(rnase1_start, rnase1_stop, rnase2_start, rnase2_stop, rnase3_start, rnase3_stop):

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
            rnase1_start = rnase1_start_list[-1]
            rnase1_stop = rnase1_stop_list[-1]
        else:
            available_slot = random.choice(rnase_slots)
            if available_slot == 'A':
                rnase1_start = random.randint(region1a_start, 122)
                rnase1_stop = rnase1_start + 10
            if available_slot == 'B':
                rnase1_start = random.randint(region1b_start, 133)
                rnase1_stop = rnase1_start + 10
            if available_slot == 'C':
                rnase1_start = random.randint(region1c_start, 149)
                rnase1_stop = rnase1_start + 10
        rnase1_start_list.append(rnase1_start)
        rnase1_stop_list.append(rnase1_stop)


    if chosen_rnase == "rnase2":
        #Adds rnase after second promoter
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
            rnase2_start = rnase2_start_list[-1]
            rnase2_stop = rnase2_stop_list[-1]
        else:
            available_slot = random.choice(rnase_slots)
            if available_slot == 'A':
                rnase2_start = random.randint(region2a_start, 281)
                rnase2_stop = rnase1_start + 10
            if available_slot == 'B':
                rnase2_start = random.randint(region2b_start, 292)
                rnase2_stop = rnase1_start + 10
            if available_slot == 'C':
                rnase2_start = random.randint(region2c_start, 308)
                rnase2_stop = rnase1_start + 10
        rnase2_start_list.append(rnase2_start)
        rnase2_stop_list.append(rnase2_stop)

    if chosen_rnase == "rnase3":
        #Adds rnase after third promoter
        rnase3_start = random.randint(region_start, (region_stop-10))
        rnase3_stop = rnase3_start + 10
        rnase3_start_list.append(rnase3_start)
        rnase3_stop_list.append(rnase3_stop)

def remove_rnase(rnase1_start, rnase1_stop, rnase2_start, rnase2_stop, rnase3_start, rnase3_stop):

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

def add_terminator(term1_start, term1_stop, term2_start, term2_stop, term3_start, term3_stop):

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
    terminator_possibilities = ["terminator1", "terminator2", "terminator3"]
    chosen_terminator = random.choice(terminator_possibilities)
    terminator_slots = ['A', 'B', 'C']

    if chosen_terminator == "terminator1":
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


    if chosen_terminator == "terminator2":
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

    if chosen_terminator == "terminator3":
        #Adds terminator after third gene
        term3_start = 449
        term3_stop = 450
        term3_start_list.append(term3_start)
        term3_stop_list.append(term3_stop)

def remove_terminator(term1_start, term1_stop, term2_start, term2_stop, term3_start, term3_stop):

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
                         rbs_start=(gene1_start-15), rbs_stop=gene1_start, rbs_strength=1e7)

        plasmid.add_gene(name="proteinY", start=gene2_start, stop=gene2_stop,
                         rbs_start=(gene2_start-15), rbs_stop=gene2_start, rbs_strength=1e7)
        plasmid.add_gene(name="proteinZ", start=gene3_start, stop=gene3_stop,
                         rbs_start=(gene3_start-15), rbs_stop=gene3_start, rbs_strength=1e7)
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
                         rbs_start=(gene1_start-15), rbs_stop=gene1_start, rbs_strength=1e7)

        plasmid.add_gene(name="proteinY", start=gene2_start, stop=gene2_stop,
                         rbs_start=(gene2_start-15), rbs_stop=gene2_start, rbs_strength=1e7)
        plasmid.add_gene(name="proteinZ", start=gene3_start, stop=gene3_stop,
                         rbs_start=(gene3_start-15), rbs_stop=gene3_start, rbs_strength=1e7)
        sim.register_genome(plasmid)
        sim.simulate(time_limit=240, time_step=1,
                     output = "best_three_genes_replicated.tsv")


if __name__ == '__main__':
    main()
