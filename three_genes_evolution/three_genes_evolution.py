import pinetree as pt
import pandas
import random
import numpy as np
import string
import filecmp

def main():

    global poly_list
    global all_poly_list
    global sos_list
    global all_sos_list
    global term_list
    global all_term_list
    sos_list = [60000]
    poly_list = []
    all_poly_list = []
    all_sos_list = []
    term_list = []
    all_term_list = []
    f_old = 1.0
    mu = 0.0
    sigma = 1.0
    poly_sigma = 1e10
    term_sigma = 1.0
    Ne = 10
    i = 0
    gen = 0

    #Taking in target tsv file
    df = pandas.read_table("three_genes_test_file.tsv", delim_whitespace=True, header=0)
    df = edit_target_file(df)
    aspect_probability = random.uniform(0.0, 1.0)

    #Evolution program
    new_pol_strength = random.randint(1e10, 5e10)
    new_term_efficiency = random.randint(1, 10)
    poly_list.append(new_pol_strength)
    all_poly_list.append(new_pol_strength)
    term_list.append(new_term_efficiency)
    all_term_list.append(new_term_efficiency)
    while i < 100:
        aspect_probability = random.uniform(0.0, 1.0)
        #Evolution program - auto changing polymerase strength
        #if aspect_probability > 0.5:
        eps = np.random.normal(mu, sigma)
        poly_eps = np.random.normal(mu, poly_sigma)
        f_new = f_old * (1.0 + eps)

        if f_new > f_old:
            #Accepting mutation
            three_genome.recreated_genome(new_pol_strength)
            #Taking in new file and removing unnecessary rows and columns
            nf = pandas.read_table("three_genes_replicated.tsv", delim_whitespace=True, header=0)
            nf = edit_new_file(nf)
            sos = sum_of_squares(df, nf)
            all_sos_list.append(sos)
            if sos == 0:
                break
            #Accepts mutation only if sum of squares value decreases
            if sos > sos_list[-1]:
                three_genome.recreated_genome(poly_list[-1])
                sos = sum_of_squares(df, nf)
            if sos <= sos_list[-1]:
                poly_list.append(new_pol_strength)
                sos_list.append(sos)
            gen+=1
        else:
            #Calculate fitness of new mutation
            probability = calc_fitness(f_new, f_old, Ne)
            f_old = probability
            if probability > random.random():
                #Accepting mutation
                three_genome.recreated_genome(new_pol_strength)
                #Taking in new file and removing unnecessary rows and columns
                nf = pandas.read_table("three_genes_replicated.tsv", delim_whitespace=True, header=0)
                nf = edit_new_file(nf)
                sos = sum_of_squares(df, nf)
                all_sos_list.append(sos)
                if sos == 0:
                    break
                #Accepts mutation only if sum of squares value decreases
                if sos > sos_list[-1]:
                    three_genome.recreated_genome(poly_list[-1])
                    sos = sum_of_squares(df, nf)
                if sos <= sos_list[-1]:
                    poly_list.append(new_pol_strength)
                    sos_list.append(sos)
                gen+=1

        #Determines new polymerase strength to reduce sum of squares value
        new_pol_strength = poly_list[-1] + poly_eps
        #if sos_list[-1] < 3000:
            #poly_sigma = 1e7
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

        #Evolution program - auto changing terminator efficiency rate
        '''if aspect_probability < 0.5:
            eps = np.random.normal(mu, sigma)
            term_eps = np.random.normal(mu, term_sigma)
            f_new = f_old * (1.0 + eps)

            if f_new > f_old:
                #Accepting mutation
                three_genome.recreated_genome(new_pol_strength, new_term_efficiency)
                #Taking in new file and removing unnecessary rows and columns
                nf = pandas.read_table("three_genes_replicated.tsv", delim_whitespace=True, header=0)
                nf = edit_new_file(nf)
                sos = sum_of_squares(df, nf)
                all_sos_list.append(sos)
                if sos == 0:
                    break
                #Accepts mutation only if sum of squares value decreases
                if sos > sos_list[-1]:
                    three_genome.recreated_genome(poly_list[-1], term_list[-1])
                    sos = sum_of_squares(df, nf)
                if sos <= sos_list[-1]:
                    term_list.append(new_term_efficiency)
                    sos_list.append(sos)
                gen+=1
            else:
                #Calculate fitness of new mutation
                probability = calc_fitness(f_new, f_old, Ne)
                f_old = probability
                if probability > random.random():
                    #Accepting mutation
                    three_genome.recreated_genome(new_pol_strengt, new_term_efficiency)
                    #Taking in new file and removing unnecessary rows and columns
                    nf = pandas.read_table("three_genes_replicated.tsv", delim_whitespace=True, header=0)
                    nf = edit_new_file(nf)
                    sos = sum_of_squares(df, nf)
                    all_sos_list.append(sos)
                    if sos == 0:
                        break
                    #Accepts mutation only if sum of squares value decreases
                    if sos > sos_list[-1]:
                        three_genome.recreated_genome(poly_list[-1], term_list[-1])
                        sos = sum_of_squares(df, nf)
                    if sos <= sos_list[-1]:
                        term_list.append(new_term_efficiency)
                        sos_list.append(sos)
                    gen+=1

            #Determines new terminator efficiency value to reduce sum of squares value
            new_term_efficiency = term_list[-1] + term_eps
            while term_list < 0:
                term_eps = np.random.normal(mu, term_sigma)
                new_term_efficiency = term_list[-1] + term_eps
                if new_term_efficiency in all_term_list:
                    term_eps = np.random.normal(mu, term_sigma)
                    new_term_efficiency = term_list[-1] + term_eps
            if new_term_efficiency in all_term_list:
                term_eps = np.random.normal(mu, term_sigma)
                new_term_efficiency = term_list[-1] + term_eps
            all_term_list.append(new_term_efficiency)'''

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
    print("Final Polymerase Rate: ", poly_list[-1])
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

    def recreated_genome(pol_strength):

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
        plasmid.add_terminator(name="t1", start=449, stop=450,
                               efficiency={"rnapol": 1.0})

        plasmid.add_gene(name="proteinX", start=26, stop=148,
                         rbs_start=(26 - 15), rbs_stop=26, rbs_strength=1e7)

        plasmid.add_gene(name="proteinY", start=26 + 150, stop=148 + 150,
                         rbs_start=(26 - 15 + 150), rbs_stop=26 + 150, rbs_strength=1e7)

        plasmid.add_gene(name="proteinZ", start=165 + 150, stop=298 + 150,
                         rbs_start=(165 - 15 + 150), rbs_stop=165 + 150, rbs_strength=1e7)
        sim.register_genome(plasmid)
        sim.simulate(time_limit=240, time_step=1,
                     output = "three_genes_replicated.tsv")




if __name__ == '__main__':
    main()
