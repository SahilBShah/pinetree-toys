import pinetree as pt
import pandas
import random
import numpy as np
import string
import filecmp

#Need to add in the following functions: remove promoters, remove rnase cleavage sites, change numeric values according to biological range, add in terminators after genes

def main():

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
    new_gene1_start = 26
    new_gene1_stop = 148
    new_gene2_start = 178
    new_gene2_stop = 300
    new_gene3_start = 315
    new_gene3_stop = 449
    f_old = 1.0
    mu = 0.0
    sigma = 1.0
    poly_sigma = 1e10
    term_sigma = 1.0
    gene_sigma = 1.0
    prom_sigma = 1.0
    Ne = 10
    i = 0
    gen = 0
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

    #Taking in target tsv file
    df = pandas.read_table(input("Enter tsv file name: "), delim_whitespace=True, header=0)
    df = edit_target_file(df)
    aspect_probability = random.uniform(0.0, 1.0)

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
    new_gene3_start = random.randint(new_gene2_stop+15, 449)
    new_gene3_stop = random.randint(new_gene3_start+1, 450)
    gene3_stop_list.append(new_gene3_stop)
    while new_gene3_start > new_gene3_stop:
        new_gene3_start = random.randint(new_gene2_stop+15, gene3_stop_list[-1])
        new_gene3_stop = random.randint(new_gene3_start+1, 450)
    new_gene2_start = random.randint(new_gene1_stop+15, 449)
    new_gene2_stop = random.randint(new_gene2_start+1, new_gene3_stop-15)
    gene2_stop_list.append(new_gene2_stop)
    while new_gene2_start > new_gene2_stop:
        new_gene2_start = random.randint(new_gene1_stop+15, gene2_stop_list[-1])
        new_gene2_stop = random.randint(new_gene2_start+1, new_gene3_stop-15)
    gene1_start_list.append(new_gene1_start)
    gene2_start_list.append(new_gene2_start)
    gene3_start_list.append(new_gene3_start)

    #Promoters that can be added in
    if new_gene2_start - (new_gene1_stop - 7) >= 25:
        new_prom1_start = random.randint(new_gene1_stop-7, new_gene2_start-25)
        new_prom1_stop = new_prom1_start + 9
    if new_gene3_start - (new_gene2_stop - 7) >= 25:
        new_prom2_start = random.randint(new_gene2_stop-7, new_gene3_start-25)
        new_prom2_stop = new_prom2_start + 9


    while i < 1000:
        aspect_probability = random.uniform(0.0, 1.0)

        if aspect_probability > 0.5:
            #Determining polymerase strength
            eps = np.random.normal(mu, sigma)
            poly_eps = np.random.normal(mu, poly_sigma)
            f_new = f_old * (1.0 + eps)

            if f_new > f_old:
                #Accepting mutation
                three_genome.recreated_genome(new_pol_strength, new_term_efficiency, new_gene1_start, new_gene1_stop,
                                                    new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop,
                                                    promoter1_start, promoter1_stop, promoter2_start, promoter2_stop,
                                                    new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop,
                                                    new_rnase3_start, new_rnase3_stop)
                #Taking in new file and removing unnecessary rows and columns
                nf = pandas.read_table("three_genes_replicated.tsv", delim_whitespace=True, header=0)
                nf = edit_new_file(nf)
                sos = sum_of_squares(df, nf)
                all_sos_list.append(sos)
                #Accepts mutation only if sum of squares value decreases
                if sos <= sos_list[-1]:
                    poly_list.append(new_pol_strength)
                    sos_list.append(sos)
                    three_genome.best_recreated_genome(new_pol_strength, new_term_efficiency, new_gene1_start, new_gene1_stop,
                                                        new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop,
                                                        promoter1_start, promoter1_stop, promoter2_start, promoter2_stop,
                                                        new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop,
                                                        new_rnase3_start, new_rnase3_stop)
                if sos_list[-1] == 0:
                    break
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
                                                        new_rnase3_start, new_rnase3_stop)
                    #Taking in new file and removing unnecessary rows and columns
                    nf = pandas.read_table("three_genes_replicated.tsv", delim_whitespace=True, header=0)
                    nf = edit_new_file(nf)
                    sos = sum_of_squares(df, nf)
                    all_sos_list.append(sos)
                    #Accepts mutation only if sum of squares value decreases
                    if sos <= sos_list[-1]:
                        poly_list.append(new_pol_strength)
                        sos_list.append(sos)
                        three_genome.best_recreated_genome(new_pol_strength, new_term_efficiency, new_gene1_start, new_gene1_stop,
                                                            new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop,
                                                            promoter1_start, promoter1_stop, promoter2_start, promoter2_stop,
                                                            new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop,
                                                            new_rnase3_start, new_rnase3_stop)
                    if sos_list[-1] == 0:
                        break
                    gen+=1

            #Determines new polymerase strength to reduce sum of squares value
            new_pol_strength = poly_list[-1] + poly_eps
            #if sos_list[-1] < 2000:
                #poly_sigma = 1e8
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


        if aspect_probability < 0.5 and aspect_probability > 0.2:
            #Determining terminator polymerase efficiency rate
            eps = np.random.normal(mu, sigma)
            term_eps = np.random.normal(mu, term_sigma)
            f_new = f_old * (1.0 + eps)

            if f_new > f_old:
                #Accepting mutation
                three_genome.recreated_genome(new_pol_strength, new_term_efficiency, new_gene1_start, new_gene1_stop,
                                                    new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop,
                                                    promoter1_start, promoter1_stop, promoter2_start, promoter2_stop,
                                                    new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop,
                                                    new_rnase3_start, new_rnase3_stop)
                #Taking in new file and removing unnecessary rows and columns
                nf = pandas.read_table("three_genes_replicated.tsv", delim_whitespace=True, header=0)
                nf = edit_new_file(nf)
                sos = sum_of_squares(df, nf)
                all_sos_list.append(sos)
                if sos == 0:
                    break
                #Accepts mutation only if sum of squares value decreases
                if sos <= sos_list[-1]:
                    term_list.append(new_term_efficiency)
                    sos_list.append(sos)
                    three_genome.best_recreated_genome(new_pol_strength, new_term_efficiency, new_gene1_start, new_gene1_stop,
                                                        new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop,
                                                        promoter1_start, promoter1_stop, promoter2_start, promoter2_stop,
                                                        new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop,
                                                        new_rnase3_start, new_rnase3_stop)
                if sos_list[-1] == 0:
                    break
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
                                                        new_rnase3_start, new_rnase3_stop)
                    #Taking in new file and removing unnecessary rows and columns
                    nf = pandas.read_table("three_genes_replicated.tsv", delim_whitespace=True, header=0)
                    nf = edit_new_file(nf)
                    sos = sum_of_squares(df, nf)
                    all_sos_list.append(sos)
                    if sos == 0:
                        break
                    #Accepts mutation only if sum of squares value decreases
                    if sos <= sos_list[-1]:
                        term_list.append(new_term_efficiency)
                        sos_list.append(sos)
                        three_genome.best_recreated_genome(new_pol_strength, new_term_efficiency, new_gene1_start, new_gene1_stop,
                                                            new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop,
                                                            promoter1_start, promoter1_stop, promoter2_start, promoter2_stop,
                                                            new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop,
                                                            new_rnase3_start, new_rnase3_stop)
                    if sos_list[-1] == 0:
                        break
                    gen+=1

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


        if aspect_probability < 0.1:
            #Changing gene lengths
            eps = np.random.normal(mu, sigma)
            gene_eps = np.random.normal(mu, gene_sigma)
            f_new = f_old * (1.0 + eps)

            if f_new > f_old:
                #Accepting mutation
                three_genome.recreated_genome(new_pol_strength, new_term_efficiency, new_gene1_start, new_gene1_stop,
                                                    new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop,
                                                    promoter1_start, promoter1_stop, promoter2_start, promoter2_stop,
                                                    new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop,
                                                    new_rnase3_start, new_rnase3_stop)
                #Taking in new file and removing unnecessary rows and columns
                nf = pandas.read_table("three_genes_replicated.tsv", delim_whitespace=True, header=0)
                nf = edit_new_file(nf)
                sos = sum_of_squares(df, nf)
                all_sos_list.append(sos)
                #Accepts mutation only if sum of squares value decreases
                if sos <= sos_list[-1]:
                    sos_list.append(sos)
                    three_genome.best_recreated_genome(new_pol_strength, new_term_efficiency, new_gene1_start, new_gene1_stop,
                                                        new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop,
                                                        promoter1_start, promoter1_stop, promoter2_start, promoter2_stop,
                                                        new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop,
                                                        new_rnase3_start, new_rnase3_stop)
                if sos_list[-1] == 0:
                    break
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
                                                        new_rnase3_start, new_rnase3_stop)
                    #Taking in new file and removing unnecessary rows and columns
                    nf = pandas.read_table("three_genes_replicated.tsv", delim_whitespace=True, header=0)
                    nf = edit_new_file(nf)
                    sos = sum_of_squares(df, nf)
                    all_sos_list.append(sos)
                    #Accepts mutation only if sum of squares value decreases
                    if sos <= sos_list[-1]:
                        sos_list.append(sos)
                        three_genome.best_recreated_genome(new_pol_strength, new_term_efficiency, new_gene1_start, new_gene1_stop,
                                                            new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop,
                                                            promoter1_start, promoter1_stop, promoter2_start, promoter2_stop,
                                                            new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop,
                                                            new_rnase3_start, new_rnase3_stop)
                    if sos_list[-1] == 0:
                        break
                    gen+=1

            #Alter gene sizes
            new_gene1_start += round(gene_eps)
            new_gene1_stop += round(gene_eps)
            new_gene2_start += round(gene_eps)
            new_gene2_stop += round(gene_eps)
            new_gene3_start += round(gene_eps)
            new_gene3_stop += round(gene_eps)

            #If gene size is larger or smaller than genome, redo
            while new_gene3_stop > 450 or new_gene1_start < 26:
                gene_eps = np.random.normal(mu, gene_sigma)
                new_gene1_start = gene1_start_list[-1] + round(gene_eps)
                new_gene1_stop = gene1_stop_list[-1] + round(gene_eps)
                new_gene2_start = gene2_start_list[-1] + round(gene_eps)
                new_gene2_stop = gene2_stop_list[-1] + round(gene_eps)
                new_gene3_start = gene3_start_list[-1] + round(gene_eps)
                new_gene3_stop = gene3_stop_list[-1] + round(gene_eps)

            gene1_start_list.append(new_gene1_start)
            gene2_start_list.append(new_gene2_start)
            gene3_start_list.append(new_gene3_start)
            gene1_stop_list.append(new_gene1_stop)
            gene2_stop_list.append(new_gene2_stop)
            gene3_stop_list.append(new_gene3_stop)


        if aspect_probability < 0.2 and aspect_probability > 0.15:
            #Adding second promoter
            if new_prom1_start > 0:
                promoter1_start = new_prom1_start
                promoter1_stop = new_prom1_stop
            eps = np.random.normal(mu, sigma)
            prom1_eps = np.random.normal(mu, prom_sigma)
            f_new = f_old * (1.0 + eps)

            if f_new > f_old:
                #Accepting mutation
                three_genome.recreated_genome(new_pol_strength, new_term_efficiency, new_gene1_start, new_gene1_stop,
                                                    new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop,
                                                    promoter1_start, promoter1_stop, promoter2_start, promoter2_stop,
                                                    new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop,
                                                    new_rnase3_start, new_rnase3_stop)
                #Taking in new file and removing unnecessary rows and columns
                nf = pandas.read_table("three_genes_replicated.tsv", delim_whitespace=True, header=0)
                nf = edit_new_file(nf)
                sos = sum_of_squares(df, nf)
                all_sos_list.append(sos)
                #Accepts mutation only if sum of squares value decreases
                if sos <= sos_list[-1]:
                    sos_list.append(sos)
                    three_genome.best_recreated_genome(new_pol_strength, new_term_efficiency, new_gene1_start, new_gene1_stop,
                                                        new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop,
                                                        promoter1_start, promoter1_stop, promoter2_start, promoter2_stop,
                                                        new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop,
                                                        new_rnase3_start, new_rnase3_stop)
                if sos_list[-1] == 0:
                    break
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
                                                        new_rnase3_start, new_rnase3_stop)
                    #Taking in new file and removing unnecessary rows and columns
                    nf = pandas.read_table("three_genes_replicated.tsv", delim_whitespace=True, header=0)
                    nf = edit_new_file(nf)
                    sos = sum_of_squares(df, nf)
                    all_sos_list.append(sos)
                    #Accepts mutation only if sum of squares value decreases
                    if sos <= sos_list[-1]:
                        sos_list.append(sos)
                        three_genome.best_recreated_genome(new_pol_strength, new_term_efficiency, new_gene1_start, new_gene1_stop,
                                                            new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop,
                                                            promoter1_start, promoter1_stop, promoter2_start, promoter2_stop,
                                                            new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop,
                                                            new_rnase3_start, new_rnase3_stop)
                    if sos_list[-1] == 0:
                        break
                    gen+=1

            new_prom1_start += round(prom1_eps)
            new_prom1_stop += round(prom1_eps)
            while new_prom1_start < new_gene2_start - 15 or new_prom1_start < new_gene1_stop - 7:
                prom1_eps = np.random.normal(mu, prom_sigma)
                new_prom1_start += round(prom1_eps)
                new_prom1_stop += round(prom1_eps)
            promoter1_start = new_prom1_start
            promoter1_stop = new_prom1_stop

        if aspect_probability < 0.15:
            #Adding third promoter
            if new_prom2_start > 0:
                promoter2_start = new_prom2_start
                promoter2_stop = new_prom2_stop
            eps = np.random.normal(mu, sigma)
            prom2_eps = np.random.normal(mu, prom_sigma)
            f_new = f_old * (1.0 + eps)

            if f_new > f_old:
                #Accepting mutation
                three_genome.recreated_genome(new_pol_strength, new_term_efficiency, new_gene1_start, new_gene1_stop,
                                                    new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop,
                                                    promoter1_start, promoter1_stop, promoter2_start, promoter2_stop,
                                                    new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop,
                                                    new_rnase3_start, new_rnase3_stop)
                #Taking in new file and removing unnecessary rows and columns
                nf = pandas.read_table("three_genes_replicated.tsv", delim_whitespace=True, header=0)
                nf = edit_new_file(nf)
                sos = sum_of_squares(df, nf)
                all_sos_list.append(sos)
                #Accepts mutation only if sum of squares value decreases
                if sos <= sos_list[-1]:
                    sos_list.append(sos)
                    three_genome.best_recreated_genome(new_pol_strength, new_term_efficiency, new_gene1_start, new_gene1_stop,
                                                        new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop,
                                                        promoter1_start, promoter1_stop, promoter2_start, promoter2_stop,
                                                        new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop,
                                                        new_rnase3_start, new_rnase3_stop)
                if sos_list[-1] == 0:
                    break
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
                                                        new_rnase3_start, new_rnase3_stop)
                    #Taking in new file and removing unnecessary rows and columns
                    nf = pandas.read_table("three_genes_replicated.tsv", delim_whitespace=True, header=0)
                    nf = edit_new_file(nf)
                    sos = sum_of_squares(df, nf)
                    all_sos_list.append(sos)
                    #Accepts mutation only if sum of squares value decreases
                    if sos <= sos_list[-1]:
                        sos_list.append(sos)
                        three_genome.best_recreated_genome(new_pol_strength, new_term_efficiency, new_gene1_start, new_gene1_stop,
                                                            new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop,
                                                            promoter1_start, promoter1_stop, promoter2_start, promoter2_stop,
                                                            new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop,
                                                            new_rnase3_start, new_rnase3_stop)
                    if sos_list[-1] == 0:
                        break
                    gen+=1

            new_prom2_start += round(prom2_eps)
            new_prom2_stop += round(prom2_eps)
            while new_prom2_start < new_gene3_start - 15 or new_prom2_start < new_gene2_stop - 7:
                prom2_eps = np.random.normal(mu, prom_sigma)
                new_prom2_start += round(prom2_eps)
                new_prom2_stop += round(prom2_eps)
            promoter2_start = new_prom2_start
            promoter2_stop = new_prom2_stop

        if aspect_probability > 0.8:
            #Rnase after initial promoter
            eps = np.random.normal(mu, sigma)
            f_new = f_old * (1.0 + eps)
            new_rnase1_start = 10
            new_rnase1_stop = 20
            if f_new > f_old:
                #Accepting mutation
                three_genome.recreated_genome(new_pol_strength, new_term_efficiency, new_gene1_start, new_gene1_stop,
                                                    new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop,
                                                    promoter1_start, promoter1_stop, promoter2_start, promoter2_stop,
                                                    new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop,
                                                    new_rnase3_start, new_rnase3_stop)
                #Taking in new file and removing unnecessary rows and columns
                nf = pandas.read_table("three_genes_replicated.tsv", delim_whitespace=True, header=0)
                nf = edit_new_file(nf)
                sos = sum_of_squares(df, nf)
                all_sos_list.append(sos)
                #Accepts mutation only if sum of squares value decreases
                if sos <= sos_list[-1]:
                    sos_list.append(sos)
                    three_genome.best_recreated_genome(new_pol_strength, new_term_efficiency, new_gene1_start, new_gene1_stop,
                                                        new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop,
                                                        promoter1_start, promoter1_stop, promoter2_start, promoter2_stop,
                                                        new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop,
                                                        new_rnase3_start, new_rnase3_stop)
                if sos_list[-1] == 0:
                    break
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
                                                        new_rnase3_start, new_rnase3_stop)
                    #Taking in new file and removing unnecessary rows and columns
                    nf = pandas.read_table("three_genes_replicated.tsv", delim_whitespace=True, header=0)
                    nf = edit_new_file(nf)
                    sos = sum_of_squares(df, nf)
                    all_sos_list.append(sos)
                    #Accepts mutation only if sum of squares value decreases
                    if sos <= sos_list[-1]:
                        sos_list.append(sos)
                        three_genome.best_recreated_genome(new_pol_strength, new_term_efficiency, new_gene1_start, new_gene1_stop,
                                                            new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop,
                                                            promoter1_start, promoter1_stop, promoter2_start, promoter2_stop,
                                                            new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop,
                                                            new_rnase3_start, new_rnase3_stop)
                    if sos_list[-1] == 0:
                        break
                    gen+=1

            if aspect_probability > 0.85:
                #Rnase after second promoter
                eps = np.random.normal(mu, sigma)
                f_new = f_old * (1.0 + eps)
                if promoter1_start > 0:
                    new_rnase1_start = promoter1_stop
                    new_rnase1_stop = promoter1_stop + 10
                if promoter1_stop > 0:
                    if f_new > f_old:
                        #Accepting mutation
                        three_genome.recreated_genome(new_pol_strength, new_term_efficiency, new_gene1_start, new_gene1_stop,
                                                            new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop,
                                                            promoter1_start, promoter1_stop, promoter2_start, promoter2_stop,
                                                            new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop,
                                                            new_rnase3_start, new_rnase3_stop)
                        #Taking in new file and removing unnecessary rows and columns
                        nf = pandas.read_table("three_genes_replicated.tsv", delim_whitespace=True, header=0)
                        nf = edit_new_file(nf)
                        sos = sum_of_squares(df, nf)
                        all_sos_list.append(sos)
                        #Accepts mutation only if sum of squares value decreases
                        if sos <= sos_list[-1]:
                            sos_list.append(sos)
                            three_genome.best_recreated_genome(new_pol_strength, new_term_efficiency, new_gene1_start, new_gene1_stop,
                                                                new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop,
                                                                promoter1_start, promoter1_stop, promoter2_start, promoter2_stop,
                                                                new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop,
                                                                new_rnase3_start, new_rnase3_stop)
                        if sos_list[-1] == 0:
                            break
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
                                                                new_rnase3_start, new_rnase3_stop)
                            #Taking in new file and removing unnecessary rows and columns
                            nf = pandas.read_table("three_genes_replicated.tsv", delim_whitespace=True, header=0)
                            nf = edit_new_file(nf)
                            sos = sum_of_squares(df, nf)
                            all_sos_list.append(sos)
                            #Accepts mutation only if sum of squares value decreases
                            if sos <= sos_list[-1]:
                                sos_list.append(sos)
                                three_genome.best_recreated_genome(new_pol_strength, new_term_efficiency, new_gene1_start, new_gene1_stop,
                                                                    new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop,
                                                                    promoter1_start, promoter1_stop, promoter2_start, promoter2_stop,
                                                                    new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop,
                                                                    new_rnase3_start, new_rnase3_stop)
                            if sos_list[-1] == 0:
                                break
                            gen+=1

            if aspect_probability > 0.90:
                #Rnase after third promoter
                eps = np.random.normal(mu, sigma)
                f_new = f_old * (1.0 + eps)
                if promoter2_start > 0:
                    new_rnase2_start = promoter2_stop
                    new_rnase2_stop = promoter2_stop + 10
                if promoter2_stop > 0:
                    if f_new > f_old:
                        #Accepting mutation
                        three_genome.recreated_genome(new_pol_strength, new_term_efficiency, new_gene1_start, new_gene1_stop,
                                                            new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop,
                                                            promoter1_start, promoter1_stop, promoter2_start, promoter2_stop,
                                                            new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop,
                                                            new_rnase3_start, new_rnase3_stop)
                        #Taking in new file and removing unnecessary rows and columns
                        nf = pandas.read_table("three_genes_replicated.tsv", delim_whitespace=True, header=0)
                        nf = edit_new_file(nf)
                        sos = sum_of_squares(df, nf)
                        all_sos_list.append(sos)
                        #Accepts mutation only if sum of squares value decreases
                        if sos <= sos_list[-1]:
                            sos_list.append(sos)
                            three_genome.best_recreated_genome(new_pol_strength, new_term_efficiency, new_gene1_start, new_gene1_stop,
                                                                new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop,
                                                                promoter1_start, promoter1_stop, promoter2_start, promoter2_stop,
                                                                new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop,
                                                                new_rnase3_start, new_rnase3_stop)
                        if sos_list[-1] == 0:
                            break
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
                                                                new_rnase3_start, new_rnase3_stop)
                            #Taking in new file and removing unnecessary rows and columns
                            nf = pandas.read_table("three_genes_replicated.tsv", delim_whitespace=True, header=0)
                            nf = edit_new_file(nf)
                            sos = sum_of_squares(df, nf)
                            all_sos_list.append(sos)
                            #Accepts mutation only if sum of squares value decreases
                            if sos <= sos_list[-1]:
                                sos_list.append(sos)
                                three_genome.best_recreated_genome(new_pol_strength, new_term_efficiency, new_gene1_start, new_gene1_stop,
                                                                    new_gene2_start, new_gene2_stop, new_gene3_start, new_gene3_stop,
                                                                    promoter1_start, promoter1_stop, promoter2_start, promoter2_stop,
                                                                    new_rnase1_start, new_rnase1_stop, new_rnase2_start, new_rnase2_stop,
                                                                    new_rnase3_start, new_rnase3_stop)
                            if sos_list[-1] == 0:
                                break
                            gen+=1


        i+=1


    #Export tsv file of sum of squares and polymerase rate data
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
    print("Final Polymerase Rate: ", poly_list[-1])
    print("Final Terminator Efficiency Rate: ", term_list[-1])
    print("Gene 1 starts: ", new_gene1_start)
    print("Gene 1 stops: ", new_gene1_stop)
    print("Gene 1 Ribosome Binding Site starts: ", new_gene1_start-15)
    print("Gene 1 Ribosome Binding Site stops: ", new_gene1_start)
    print("Gene 2 starts: ", new_gene2_start)
    print("Gene 2 stops: ", new_gene2_stop)
    print("Gene 2 Ribosome Binding Site starts: ", new_gene2_start-15)
    print("Gene 2 Ribosome Binding Site stops: ", new_gene2_start)
    print("Gene 3 starts: ", new_gene3_start)
    print("Gene 3 stops: ", new_gene3_stop)
    print("Gene 3 Ribosome Binding Site starts: ", new_gene3_start-15)
    print("Gene 3 Ribosome Binding Site stops: ", new_gene3_start)
    print("Added Promoter 1 starts at: ", promoter1_start)
    print("Added Promoter 1 stops at: ", promoter1_stop)
    print("Added Promoter 2 starts at: ", promoter2_start)
    print("Added Promoter 2 stops at: ", promoter2_stop)
    print("Rnase site 1 starts at: ", new_rnase1_start)
    print("Rnase site 1 stops at: ", new_rnase1_stop)
    print("Rnase site 2 starts at: ", new_rnase2_start)
    print("Rnase site 2 stops at: ", new_rnase2_stop)
    print("Rnase site 3 starts at: ", new_rnase3_start)
    print("Rnase site 3 stops at: ", new_rnase3_stop)
    print("Generations = ", gen)


#Calculates the fitness of the new mutation
def calc_fitness(variant_fit, orig_fit, Ne):

    xi = orig_fit
    xj = variant_fit
    N = Ne

    if xj == xi:
        return(1.0 / float(N) )
    if xj > xi:
        return(1.0)
    else:
        variant_fit = (xj / xi) ** (2 * N - 2)
        return variant_fit

    try:
        p =((1-pow((xi/xj), 2)) /(1-pow((xi/xj), (2 * float(N) )) ) )
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

#Class containing genome for simulation with new polymerase strength
class three_genome:

    def __init__(self):
        self.pol_strength = 4e10

    def recreated_genome(pol_strength, term_efficiency, gene1_start, gene1_stop, gene2_start,
                         gene2_stop, gene3_start, gene3_stop, prom1_start, prom1_stop, prom2_start,
                         prom2_stop, rnase1_start, rnase1_stop, rnase2_start, rnase2_stop, rnase3_start,
                         rnase3_stop):

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
        plasmid.add_terminator(name="t1", start=449, stop=450,
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
                              rnase3_stop):

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
        plasmid.add_terminator(name="t1", start=449, stop=450,
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
