import pinetree as pt
import pandas
import random
import numpy as np
import string
import filecmp
import math
import os

#Need to add in the following functions: change numeric values according to biological range, expand genome when adding in components

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
    gene1_list = []
    gene2_list = []
    gene3_list = []
    term1_start_list = [0]
    term1_stop_list = [0]
    term2_start_list = [0]
    term2_stop_list = [0]
    term3_start_list = [0]
    term3_stop_list = [0]
    term1_list = []
    term2_list = []
    term3_list = []
    rnase1_start_list = [0]
    rnase1_stop_list = [0]
    rnase2_start_list = [0]
    rnase2_stop_list = [0]
    rnase3_start_list = [0]
    rnase3_stop_list = [0]
    rnase1_list = []
    rnase2_list = []
    rnase3_list = []
    prom1_start_list = [0]
    prom1_stop_list = [0]
    prom2_start_list = [0]
    prom2_stop_list = [0]
    prom1_list = []
    prom2_list = []
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

    '''possibilities = ["alter polymerase 1 strength", "alter terminator efficiency", "alter gene length", "add promoter", "remove promoter", "add rnase",
                        "remove rnase", "add terminator", "remove terminator"]'''

    possibilities = ["alter polymerase 1 strength", "add promoter", "remove promoter", "add rnase",
                      "remove rnase", "add terminator", "remove terminator"]

    while i < 100:

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
            prom1_start_list = new_promoter[0]
            prom1_stop_list = new_promoter[1]
            prom2_start_list = new_promoter[2]
            prom2_stop_list = new_promoter[3]
            new_prom1_start = prom1_start_list[-1]
            new_prom1_stop = prom1_stop_list[-1]
            new_prom2_start = prom2_start_list[-1]
            new_prom2_stop = prom2_stop_list[-1]
            new_pol1_strength = new_promoter[-2]
            new_pol2_strength = new_promoter[-1]

        if mutation == "remove promoter":
            old_promoter = remove_promoter(new_prom1_start, new_prom1_stop, new_prom2_start, new_prom2_stop, prom1_start_list, prom1_stop_list,
                                            prom2_start_list, prom2_stop_list)
            prom1_start_list = old_promoter[0]
            prom1_stop_list = old_promoter[1]
            prom2_start_list = old_promoter[2]
            prom2_stop_list = old_promoter[3]
            new_prom1_start = prom1_start_list[-1]
            new_prom1_stop = prom1_stop_list[-1]
            new_prom2_start = prom2_start_list[-1]
            new_prom2_stop = prom2_stop_list[-1]

        if mutation == "add rnase":
            new_rnase = add_ranse(new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop, new_rnase3_start, new_rnase3_stop,
                                    prom1_start_list, prom1_stop_list, prom2_start_list, prom2_stop_list, term1_start_list, term1_stop_list,
                                    term2_start_list, term2_stop_list, rnase1_start_list, rnase1_stop_list, rnase2_start_list, rnase2_stop_list,
                                    rnase3_start_list, rnase3_stop_list)
            rnase1_start_list = new_rnase[0]
            rnase1_stop_list = new_rnase[1]
            rnase2_start_list = new_rnase[2]
            rnase2_stop_list = new_rnase[3]
            rnase3_start_list = new_rnase[4]
            rnase3_stop_list = new_rnase[5]
            new_rnase1_start = rnase1_start_list[-1]
            new_rnase1_stop = rnase1_stop_list[-1]
            new_rnase2_start = rnase2_start_list[-1]
            new_rnase2_stop = rnase2_stop_list[-1]
            new_rnase3_start = rnase3_start_list[-1]
            new_rnase3_stop = rnase3_stop_list[-1]

        if mutation == "remove rnase":
            old_rnase = remove_rnase(new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop, new_rnase3_start, new_rnase3_stop,
                                        rnase1_start_list, rnase1_stop_list, rnase2_start_list, rnase2_stop_list, rnase3_start_list, rnase3_stop_list)
            rnase1_start_list = old_rnase[0]
            rnase1_stop_list = old_rnase[1]
            rnase2_start_list = old_rnase[2]
            rnase2_stop_list = old_rnase[3]
            rnase3_start_list = old_rnase[4]
            rnase3_stop_list = old_rnase[5]
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
            term1_start_list = new_terminator[0]
            term1_stop_list = new_terminator[1]
            term2_start_list = new_terminator[2]
            term2_stop_list = new_terminator[3]
            term3_start_list = new_terminator[4]
            term3_stop_list = new_terminator[5]
            new_term1_start = term1_start_list[-1]
            new_term1_stop = term1_stop_list[-1]
            new_term2_start = term2_start_list[-1]
            new_term2_stop = term2_stop_list[-1]
            new_term3_start = term3_start_list[-1]
            new_term3_stop = term3_stop_list[-1]
            new_term1_efficiency = new_terminator[-3]
            new_term2_efficiency = new_terminator[-2]
            new_term3_efficiency = new_terminator[-1]

        if mutation == "remove terminator":
            old_terminator = remove_terminator(new_term1_start, new_term1_stop, new_term2_start, new_term2_stop, new_term3_start, new_term3_stop,
                                                    term1_start_list, term1_stop_list, term2_start_list, term2_stop_list, term3_start_list, term3_stop_list)
            term1_start_list = old_terminator[0]
            term1_stop_list = old_terminator[1]
            term2_start_list = old_terminator[2]
            term2_stop_list = old_terminator[3]
            term3_start_list = old_terminator[4]
            term3_stop_list = old_terminator[5]
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
                                    rnase1_list, rnase2_list, rnase3_list)
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
            gene1_list.append(new_gene1_start)
            gene1_list.append(new_gene1_stop)
            gene2_list.append(new_gene2_start)
            gene2_list.append(new_gene2_stop)
            gene3_list.append(new_gene3_start)
            gene3_list.append(new_gene3_stop)
            prom1_list.append(new_prom1_start)
            prom1_list.append(new_prom1_stop)
            prom2_list.append(new_prom2_start)
            prom2_list.append(new_prom2_stop)
            term1_list.append(new_term1_start)
            term1_list.append(new_term1_stop)
            term2_list.append(new_term2_start)
            term2_list.append(new_term2_stop)
            term3_list.append(new_term3_start)
            term3_list.append(new_term3_stop)
            rnase1_list.append(new_rnase1_start)
            rnase1_list.append(new_rnase1_stop)
            rnase2_list.append(new_rnase2_start)
            rnase2_list.append(new_rnase2_stop)
            rnase3_list.append(new_rnase3_start)
            rnase3_list.append(new_rnase3_stop)
            get_coordinates(gene1_list, gene2_list, gene3_list, prom1_list, prom2_list, term1_list, term2_list, term3_list,
                                rnase1_list, rnase2_list, rnase3_list, iteration)

        if i == 500:
            output_file = "gen_500_data.tsv"
            three_genome.recreated_genome(new_pol_strength, new_pol1_strength, new_pol2_strength, new_term1_efficiency,
                                            new_term2_efficiency, new_term3_efficiency, new_gene1_start, new_gene1_stop,
                                            new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop,
                                            new_prom1_start, new_prom1_stop, new_prom2_start, new_prom2_stop,
                                            new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop,
                                            new_rnase3_start, new_rnase3_stop, new_term1_start, new_term1_stop,
                                            new_term2_start, new_term2_stop, new_term3_start, new_term3_stop, output_file)
            #Gets genome coordinates at the 500th iteration
            iteration = "500"
            gene1_list.append(new_gene1_start)
            gene1_list.append(new_gene1_stop)
            gene2_list.append(new_gene2_start)
            gene2_list.append(new_gene2_stop)
            gene3_list.append(new_gene3_start)
            gene3_list.append(new_gene3_stop)
            prom1_list.append(new_prom1_start)
            prom1_list.append(new_prom1_stop)
            prom2_list.append(new_prom2_start)
            prom2_list.append(new_prom2_stop)
            term1_list.append(new_term1_start)
            term1_list.append(new_term1_stop)
            term2_list.append(new_term2_start)
            term2_list.append(new_term2_stop)
            term3_list.append(new_term3_start)
            term3_list.append(new_term3_stop)
            rnase1_list.append(new_rnase1_start)
            rnase1_list.append(new_rnase1_stop)
            rnase2_list.append(new_rnase2_start)
            rnase2_list.append(new_rnase2_stop)
            rnase3_list.append(new_rnase3_start)
            rnase3_list.append(new_rnase3_stop)
            get_coordinates(gene1_list, gene2_list, gene3_list, prom1_list, prom2_list, term1_list, term2_list, term3_list,
                                rnase1_list, rnase2_list, rnase3_list, iteration)

        if i == 999:
            output_file = "gen_1000_data.tsv"
            three_genome.recreated_genome(new_pol_strength, new_pol1_strength, new_pol2_strength, new_term1_efficiency,
                                            new_term2_efficiency, new_term3_efficiency, new_gene1_start, new_gene1_stop,
                                            new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop,
                                            new_prom1_start, new_prom1_stop, new_prom2_start, new_prom2_stop,
                                            new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop,
                                            new_rnase3_start, new_rnase3_stop, new_term1_start, new_term1_stop,
                                            new_term2_start, new_term2_stop, new_term3_start, new_term3_stop, output_file)
            #Gets genome coordinates at the 1000th iteration
            iteration = "999"
            gene1_list.append(new_gene1_start)
            gene1_list.append(new_gene1_stop)
            gene2_list.append(new_gene2_start)
            gene2_list.append(new_gene2_stop)
            gene3_list.append(new_gene3_start)
            gene3_list.append(new_gene3_stop)
            prom1_list.append(new_prom1_start)
            prom1_list.append(new_prom1_stop)
            prom2_list.append(new_prom2_start)
            prom2_list.append(new_prom2_stop)
            term1_list.append(new_term1_start)
            term1_list.append(new_term1_stop)
            term2_list.append(new_term2_start)
            term2_list.append(new_term2_stop)
            term3_list.append(new_term3_start)
            term3_list.append(new_term3_stop)
            rnase1_list.append(new_rnase1_start)
            rnase1_list.append(new_rnase1_stop)
            rnase2_list.append(new_rnase2_start)
            rnase2_list.append(new_rnase2_stop)
            rnase3_list.append(new_rnase3_start)
            rnase3_list.append(new_rnase3_stop)
            get_coordinates(gene1_list, gene2_list, gene3_list, prom1_list, prom2_list, term1_list, term2_list, term3_list,
                                rnase1_list, rnase2_list, rnase3_list, iteration)

        i+=1
        print("i = ", i)


    #Exported tsv files
    '''iteration+=1
    sim_number = iteration
    path = "~/pinetree-toys/three_genes_evolution/" + string(sim_number)
    if not os.path.exists(path):
        os.makedirs(path)'''
    all_sos_dataframe = pandas.DataFrame(all_sos_list, columns=["Sum_of_Squares"])
    export_csv = all_sos_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/all_sos_data.tsv", index=False)
    sos_dataframe = pandas.DataFrame(sos_list, columns=["Sum_of_Squares"])
    export_csv = sos_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/sos_data.tsv", index=False)
    all_poly_dataframe = pandas.DataFrame(all_poly_list, columns=["Polymerase_Rate"])
    export_csv = all_poly_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/all_poly_data.tsv", index=False)
    poly_dataframe = pandas.DataFrame(poly_list, columns=["Polymerase_Rate"])
    export_csv = poly_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/poly_data.tsv", index=False)
    poly1_dataframe = pandas.DataFrame(poly1_list, columns=["Polymerase_Rate"])
    export_csv = poly1_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/poly1_data.tsv", index=False)
    poly2_dataframe = pandas.DataFrame(poly2_list, columns=["Polymerase_Rate"])
    export_csv = poly2_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/poly2_data.tsv", index=False)
    all_term_dataframe = pandas.DataFrame(all_term_list, columns=["Terminator Efficiency Rate"])
    export_csv = all_term_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/all_term_data.tsv", index=False)
    term1_efficiency_dataframe = pandas.DataFrame(term1_efficiency_list, columns=["Terminator1_Efficiency_Rate"])
    export_csv = term1_efficiency_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/term1_efficiency_data.tsv", index=False)
    term2_efficiency_dataframe = pandas.DataFrame(term2_efficiency_list, columns=["Terminator2_Efficiency_Rate"])
    export_csv = term2_efficiency_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/term2_efficiency_data.tsv", index=False)
    term3_efficiency_dataframe = pandas.DataFrame(term3_efficiency_list, columns=["Terminator3_Efficiency_Rate"])
    export_csv = term3_efficiency_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/term3_efficiency_data.tsv", index=False)
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
    print("Final Polymerase 1 Rate: ", poly_list[-1])
    print("Final Polymerase 2 Rate: ", poly1_list[-1])
    print("Final Polymerase 3 Rate: ", poly2_list[-1])
    print("Final Terminator 1 Efficiency Rate: ", term1_efficiency_list[-1])
    print("Final Terminator 2 Efficiency Rate: ", term2_efficiency_list[-1])
    print("Final Terminator 3 Efficiency Rate: ", term3_efficiency_list[-1])
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
    print("Rnase site 1 starts at: ", rnase1_start_list[-1], "and stops at: ", rnase1_stop_list[-1])
    print("Rnase site 2 starts at: ", rnase2_start_list[-1], "and stops at: ", rnase2_stop_list[-1])
    print("Rnase site 3 starts at: ", rnase3_start_list[-1], "and stops at: ", rnase3_stop_list[-1])
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

def accept_mutation(df, f_old, f_new, new_pol_strength, new_pol1_strength, new_pol2_strength, new_term1_efficiency, new_term2_efficiency,
                    new_term3_efficiency, new_gene1_start, new_gene1_stop, new_gene2_start,
                    new_gene2_stop, new_gene3_start, new_gene3_stop, promoter1_start, promoter1_stop, promoter2_start,
                    promoter2_stop, new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop, new_rnase3_start,
                    new_rnase3_stop, new_term1_start, new_term1_stop, new_term2_start, new_term2_stop, new_term3_start,
                    new_term3_stop, all_sos_list, sos_list, poly_list, poly1_list, poly2_list, term1_efficiency_list, term2_efficiency_list,
                    term3_efficiency_list, gene1_list, gene2_list, gene3_list, prom1_list, prom2_list, term1_list, term2_list, term3_list
                    rnase1_list, rnase2_list, rnase3_list):

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
        nf = pandas.read_table("three_genes_replicated.tsv", delim_whitespace=True, header=0)
        nf = edit_new_file(nf)
        sos = sum_of_squares(df, nf)
        all_sos_list.append(sos)
        #Accepts mutation only if sum of squares value decreases
        if sos <= sos_list[-1]:
            poly_list.append(new_pol_strength)
            poly1_list.append(new_pol1_strength)
            poly2_list.append(new_pol2_strength)
            term1_efficiency_list.append(new_term1_efficiency)
            term2_efficiency_list.append(new_term2_efficiency)
            term3_efficiency_list.append(new_term3_efficiency)
            sos_list.append(sos)
            three_genome.best_recreated_genome(new_pol_strength, new_pol1_strength, new_pol2_strength, new_term1_efficiency,
                                                new_term2_efficiency, new_term3_efficiency, new_gene1_start, new_gene1_stop,
                                                new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop,
                                                promoter1_start, promoter1_stop, promoter2_start, promoter2_stop,
                                                new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop,
                                                new_rnase3_start, new_rnase3_stop, new_term1_start, new_term1_stop,
                                                new_term2_start, new_term2_stop, new_term3_start, new_term3_stop)
            #Gets genome coordinates
            gene1_list.append(new_gene1_start)
            gene1_list.append(new_gene1_stop)
            gene2_list.append(new_gene2_start)
            gene2_list.append(new_gene2_stop)
            gene3_list.append(new_gene3_start)
            gene3_list.append(new_gene3_stop)
            prom1_list.append(new_prom1_start)
            prom1_list.append(new_prom1_stop)
            prom2_list.append(new_prom2_start)
            prom2_list.append(new_prom2_stop)
            term1_list.append(new_term1_start)
            term1_list.append(new_term1_stop)
            term2_list.append(new_term2_start)
            term2_list.append(new_term2_stop)
            term3_list.append(new_term3_start)
            term3_list.append(new_term3_stop)
            rnase1_list.append(new_rnase1_start)
            rnase1_list.append(new_rnase1_stop)
            rnase2_list.append(new_rnase2_start)
            rnase2_list.append(new_rnase2_stop)
            rnase3_list.append(new_rnase3_start)
            rnase3_list.append(new_rnase3_stop)
            get_coordinates(gene1_list, gene2_list, gene3_list, prom1_list, prom2_list, term1_list, term2_list, term3_list,
                                rnase1_list, rnase2_list, rnase3_list, iteration)
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
            nf = pandas.read_table("three_genes_replicated.tsv", delim_whitespace=True, header=0)
            nf = edit_new_file(nf)
            sos = sum_of_squares(df, nf)
            all_sos_list.append(sos)
            #Accepts mutation only if sum of squares value decreases
            if sos <= sos_list[-1]:
                poly_list.append(new_pol_strength)
                poly1_list.append(new_pol1_strength)
                poly2_list.append(new_pol2_strength)
                term1_efficiency_list.append(new_term1_efficiency)
                term2_efficiency_list.append(new_term2_efficiency)
                term3_efficiency_list.append(new_term3_efficiency)
                sos_list.append(sos)
                three_genome.best_recreated_genome(new_pol_strength, new_pol1_strength, new_pol2_strength, new_term1_efficiency,
                                                    new_term2_efficiency, new_term3_efficiency, new_gene1_start, new_gene1_stop,
                                                    new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop,
                                                    promoter1_start, promoter1_stop, promoter2_start, promoter2_stop,
                                                    new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop,
                                                    new_rnase3_start, new_rnase3_stop, new_term1_start, new_term1_stop,
                                                    new_term2_start, new_term2_stop, new_term3_start, new_term3_stop)
                #Gets genome coordinates
                gene1_list.append(new_gene1_start)
                gene1_list.append(new_gene1_stop)
                gene2_list.append(new_gene2_start)
                gene2_list.append(new_gene2_stop)
                gene3_list.append(new_gene3_start)
                gene3_list.append(new_gene3_stop)
                prom1_list.append(new_prom1_start)
                prom1_list.append(new_prom1_stop)
                prom2_list.append(new_prom2_start)
                prom2_list.append(new_prom2_stop)
                term1_list.append(new_term1_start)
                term1_list.append(new_term1_stop)
                term2_list.append(new_term2_start)
                term2_list.append(new_term2_stop)
                term3_list.append(new_term3_start)
                term3_list.append(new_term3_stop)
                rnase1_list.append(new_rnase1_start)
                rnase1_list.append(new_rnase1_stop)
                rnase2_list.append(new_rnase2_start)
                rnase2_list.append(new_rnase2_stop)
                rnase3_list.append(new_rnase3_start)
                rnase3_list.append(new_rnase3_stop)
                get_coordinates(gene1_list, gene2_list, gene3_list, prom1_list, prom2_list, term1_list, term2_list, term3_list,
                                    rnase1_list, rnase2_list, rnase3_list, iteration)
                gen+=1
    if probability > 0:
        return probability
    else:
        return 0

def alter_poly_strength(poly_list, poly1_list, poly2_list, promoter):

    poly_sigma = 1e10

    if promoter == "promoter":
        #Determining polymerase strength
        poly_eps = np.random.normal(mu, poly_sigma)
        #Determines new polymerase strength to reduce sum of squares value
        pol_strength = poly_list[-1] + poly_eps
        while pol_strength < 0:
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
        return pol1_strength

    if promoter == "promoter2":
        #Determining polymerase strength
        poly_eps = np.random.normal(mu, poly_sigma)
        #Determines new polymerase strength to reduce sum of squares value
        pol2_strength = poly2_list[-1] + poly_eps
        while pol2_strength < 0:
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

'''def alter_gene_length(gene1_start, gene1_stop, gene2_start, gene2_stop, gene3_start, gene3_stop):

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
        gene3_stop_list.append(gene3_stop)'''


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

    prom_list = [prom1_start_list, prom1_stop_list, prom2_start_list, prom2_stop_list, pol1_strength, pol2_strength]
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
    prom_list = [prom1_start_list, prom1_stop_list, prom2_start_list, prom2_stop_list]
    return prom_list

def add_ranse(rnase1_start, rnase1_stop, rnase2_start, rnase2_stop, rnase3_start, rnase3_stop, prom1_start_list,
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
    rnase_list = [rnase1_start_list, rnase1_stop_list, rnase2_start_list, rnase2_stop_list, rnase3_start_list, rnase3_stop_list]
    return rnase_list

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
    rnase_list = [rnase1_start_list, rnase1_stop_list, rnase2_start_list, rnase2_stop_list, rnase3_start_list, rnase3_stop_list]
    return rnase_list

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

    term_list = [term1_start_list, term1_stop_list, term2_start_list, term2_stop_list, term3_start_list, term3_stop_list, term1_efficiency, term2_efficiency, term3_efficiency]
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
    term_list = [term1_start_list, term1_stop_list, term2_start_list, term2_stop_list, term3_start_list, term3_stop_list]
    return term_list

def get_coordinates(gene1_list, gene2_list, gene3_list, prom1_list, prom2_list, term1_list, term2_list, term3_list, rnase1_list, rnase2_list, rnase3_list, generation):

    #Gets coordinates of genome and outputs as tsv file
    if generation == "best":
        gene1_dataframe = pandas.DataFrame(gene1_list, columns=["Gene_Location"])
        export_csv = gene1_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/best_genome/gene1_data.tsv", index=False)
        gene2_dataframe = pandas.DataFrame(gene2_list, columns=["Gene_Location"])
        export_csv = gene2_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/best_genome/gene2_data.tsv", index=False)
        gene3_dataframe = pandas.DataFrame(gene3_list, columns=["Gene_Location"])
        export_csv = gene3_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/best_genome/gene3_data.tsv", index=False)
        prom1_dataframe = pandas.DataFrame(prom1_list, columns=["Promoter_Location"])
        export_csv = prom1_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/best_genome/promoter1_data.tsv", index=False)
        prom2_dataframe = pandas.DataFrame(prom1_list, columns=["Promoter_Location"])
        export_csv = prom2_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/best_genome/promoter2_data.tsv", index=False)
        term1_dataframe = pandas.DataFrame(term1_list, columns=["Terminator_Location"])
        export_csv = term1_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/best_genome/term1_data.tsv", index=False)
        term2_dataframe = pandas.DataFrame(term2_list, columns=["Terminator_Location"])
        export_csv = term2_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/best_genome/term2_data.tsv", index=False)
        term3_dataframe = pandas.DataFrame(term3_list, columns=["Terminator_Location"])
        export_csv = term3_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/best_genome/term3_data.tsv", index=False)
        rnase1_dataframe = pandas.DataFrame(rnase1_list, columns=["Rnase_Location"])
        export_csv = rnase1_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/best_genome/rnase1_data.tsv", index=False)
        rnase2_dataframe = pandas.DataFrame(rnase2_list, columns=["Rnase_Location"])
        export_csv = rnase2_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/best_genome/rnase2_data.tsv", index=False)
        rnase3_dataframe = pandas.DataFrame(rnase3_list, columns=["Rnase_Location"])
        export_csv = rnase3_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/best_genome/rnase3_data.tsv", index=False)

    if generation == "0":
        gene1_dataframe = pandas.DataFrame(gene1_list, columns=["Gene_Location"])
        export_csv = gene1_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen0/gene1_data.tsv", index=False)
        gene2_dataframe = pandas.DataFrame(gene2_list, columns=["Gene_Location"])
        export_csv = gene2_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen0/gene2_data.tsv", index=False)
        gene3_dataframe = pandas.DataFrame(gene3_list, columns=["Gene_Location"])
        export_csv = gene3_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen0/gene3_data.tsv", index=False)
        prom1_dataframe = pandas.DataFrame(prom1_list, columns=["Promoter_Location"])
        export_csv = prom1_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen0/promoter1_data.tsv", index=False)
        prom2_dataframe = pandas.DataFrame(prom1_list, columns=["Promoter_Location"])
        export_csv = prom2_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen0/promoter2_data.tsv", index=False)
        term1_dataframe = pandas.DataFrame(term1_list, columns=["Terminator_Location"])
        export_csv = term1_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen0/term1_data.tsv", index=False)
        term2_dataframe = pandas.DataFrame(term2_list, columns=["Terminator_Location"])
        export_csv = term2_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen0/term2_data.tsv", index=False)
        term3_dataframe = pandas.DataFrame(term3_list, columns=["Terminator_Location"])
        export_csv = term3_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen0/term3_data.tsv", index=False)
        rnase1_dataframe = pandas.DataFrame(rnase1_list, columns=["Rnase_Location"])
        export_csv = rnase1_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen0/rnase1_data.tsv", index=False)
        rnase2_dataframe = pandas.DataFrame(rnase2_list, columns=["Rnase_Location"])
        export_csv = rnase2_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen0/rnase2_data.tsv", index=False)
        rnase3_dataframe = pandas.DataFrame(rnase3_list, columns=["Rnase_Location"])
        export_csv = rnase3_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen0/rnase3_data.tsv", index=False)

    if generation == "500":
        gene1_dataframe = pandas.DataFrame(gene1_list, columns=["Gene_Location"])
        export_csv = gene1_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen500/gene1_data.tsv", index=False)
        gene2_dataframe = pandas.DataFrame(gene2_list, columns=["Gene_Location"])
        export_csv = gene2_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen500/gene2_data.tsv", index=False)
        gene3_dataframe = pandas.DataFrame(gene3_list, columns=["Gene_Location"])
        export_csv = gene3_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen500/gene3_data.tsv", index=False)
        prom1_dataframe = pandas.DataFrame(prom1_list, columns=["Promoter_Location"])
        export_csv = prom1_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen500/promoter1_data.tsv", index=False)
        prom2_dataframe = pandas.DataFrame(prom1_list, columns=["Promoter_Location"])
        export_csv = prom2_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen500/promoter2_data.tsv", index=False)
        term1_dataframe = pandas.DataFrame(term1_list, columns=["Terminator_Location"])
        export_csv = term1_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen500/term1_data.tsv", index=False)
        term2_dataframe = pandas.DataFrame(term2_list, columns=["Terminator_Location"])
        export_csv = term2_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen500/term2_data.tsv", index=False)
        term3_dataframe = pandas.DataFrame(term3_list, columns=["Terminator_Location"])
        export_csv = term3_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen500/term3_data.tsv", index=False)
        rnase1_dataframe = pandas.DataFrame(rnase1_list, columns=["Rnase_Location"])
        export_csv = rnase1_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen500/rnase1_data.tsv", index=False)
        rnase2_dataframe = pandas.DataFrame(rnase2_list, columns=["Rnase_Location"])
        export_csv = rnase2_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen500/rnase2_data.tsv", index=False)
        rnase3_dataframe = pandas.DataFrame(rnase3_list, columns=["Rnase_Location"])
        export_csv = rnase3_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen500/rnase3_data.tsv", index=False)

    if generation == "1000":
        gene1_dataframe = pandas.DataFrame(gene1_list, columns=["Gene_Location"])
        export_csv = gene1_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen1000/gene1_data.tsv", index=False)
        gene2_dataframe = pandas.DataFrame(gene2_list, columns=["Gene_Location"])
        export_csv = gene2_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen1000/gene2_data.tsv", index=False)
        gene3_dataframe = pandas.DataFrame(gene3_list, columns=["Gene_Location"])
        export_csv = gene3_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen1000/gene3_data.tsv", index=False)
        prom1_dataframe = pandas.DataFrame(prom1_list, columns=["Promoter_Location"])
        export_csv = prom1_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen1000/promoter1_data.tsv", index=False)
        prom2_dataframe = pandas.DataFrame(prom1_list, columns=["Promoter_Location"])
        export_csv = prom2_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen1000/promoter2_data.tsv", index=False)
        term1_dataframe = pandas.DataFrame(term1_list, columns=["Terminator_Location"])
        export_csv = term1_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen1000/term1_data.tsv", index=False)
        term2_dataframe = pandas.DataFrame(term2_list, columns=["Terminator_Location"])
        export_csv = term2_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen1000/term2_data.tsv", index=False)
        term3_dataframe = pandas.DataFrame(term3_list, columns=["Terminator_Location"])
        export_csv = term3_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen1000/term3_data.tsv", index=False)
        rnase1_dataframe = pandas.DataFrame(rnase1_list, columns=["Rnase_Location"])
        export_csv = rnase1_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen1000/rnase1_data.tsv", index=False)
        rnase2_dataframe = pandas.DataFrame(rnase2_list, columns=["Rnase_Location"])
        export_csv = rnase2_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen1000/rnase2_data.tsv", index=False)
        rnase3_dataframe = pandas.DataFrame(rnase3_list, columns=["Rnase_Location"])
        export_csv = rnase3_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/gen1000/rnase3_data.tsv", index=False)

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
                                term1_stop, term2_start, term2_stop, term3_start, term3_stop):

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
                     output = "best_three_genes_replicated.tsv")


if __name__ == '__main__':
    main()
