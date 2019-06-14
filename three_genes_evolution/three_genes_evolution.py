import pinetree as pt
import pandas
import random
import numpy as np
import string
import filecmp

def main():

    global poly_list
    global sos_list
    sos_list = [200000]
    poly_list = []
    f_old = 1.0
    mu = 0.0
    sigma = 1.0
    Ne = 10
    i = 0
    gen = 0

    #Taking in target tsv file
    df = pandas.read_table("three_genes_rnase.tsv", delim_whitespace=True, header=0)
    df = edit_target_file(df)

    #Creating starting three genes sequence
    sim = pt.Model(cell_volume=8e-16)
    sim.seed(34)
    sim.add_polymerase(name="rnapol", copy_number=4, speed=40, footprint=10)
    sim.add_ribosome(copy_number=100, speed=30, footprint=10)

    plasmid = pt.Genome(name="plasmid", length=450,
                        transcript_degradation_rate=1e-2,
                        transcript_degradation_rate_ext=1e-2,
                        rnase_speed=20,
                        rnase_footprint=10)

    # plasmid = pt.Genome(name="plasmid", length=300)

    plasmid.add_promoter(name="p1", start=1, stop=10,
                         interactions={"rnapol": 4e10})
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
                 output = "three_genes_rnase2.tsv")

    #Evolution program - auto changing polymerase strength
    new_pol_strength = random.randint(1e10, 5e10)
    poly_list.append(new_pol_strength)
    while i < 1000:
        eps = np.random.normal(mu, sigma)
        f_new = f_old * (1.0 + eps)
        #new_pol_strength = random.randint(1e10, 5e10)

        if f_new > f_old:
            #Accepting mutation
            three_genome.recreated_genome(new_pol_strength)
            #Taking in new file and removing unnecessary rows and columns
            nf = pandas.read_table("three_genes_rnase2.tsv", delim_whitespace=True, header=0)
            nf = edit_new_file(nf)
            sos = sum_of_squares(df, nf)
            #Accepts mutation only if sum of squares value decreases
            if sos > sos_list[-1]:
                three_genome.recreated_genome(poly_list[-1])
                sos = sum_of_squares(df, nf)
            if sos < sos_list[-1]:
                poly_list.append(new_pol_strength)
                sos_list.append(sos)
            print("Poly List = ", poly_list)
            print("SOS List = ", sos_list)
            gen+=1
        else:
            #Calculate fitness of new mutation
            probability = calc_fitness(f_new, f_old, Ne)
            f_old = probability
            if probability > random.random():
                #Accepting mutation
                three_genome.recreated_genome(new_pol_strength)
                #Taking in new file and removing unnecessary rows and columns
                nf = pandas.read_table("three_genes_rnase2.tsv", delim_whitespace=True, header=0)
                nf = edit_new_file(nf)
                sos = sum_of_squares(df, nf)
                #Accepts mutation only if sum of squares value decreases
                if sos > sos_list[-1]:
                    three_genome.recreated_genome(poly_list[-1])
                    sos = sum_of_squares(df, nf)
                if sos < sos_list[-1]:
                    poly_list.append(new_pol_strength)
                    sos_list.append(sos)
                    new_pol_strength = random.randint(poly_list[-1], poly_list[-2])
                print("Poly List = ", poly_list)
                print("SOS List = ", sos_list)
                gen+=1

        #Determines new polymerase strength to reduce sum of squares value
        new_pol_strength = random.randint(1e10, poly_list[-1])
        i+=1
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
        p =((1-pow( (xi/xj),2)) /(1-pow( (xi/xj), (2 * float(N) )) ) )
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
                     output = "three_genes_rnase2.tsv")




if __name__ == '__main__':
    main()
