import pinetree as pt
import pandas
import random
import numpy as np
import string

#Taking in target tsv file
#df = pandas.read_csv(input("Enter .tsv file name: "), sep="/t", engine="python")
df = pandas.read_csv("three_genes_rnase.tsv", sep="/t", engine="python")

def main():

    f_old = 1.0
    mu = 0
    sigma = 1
    Ne = 10
    i = 0

    #Creating starting three genes sequence
    '''sim = pt.Model(cell_volume=8e-16)
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
                 output = "three_genes_rnase2.tsv")'''

    #Evolution program - auto changing polymerase strength
    prom = 0
    while i < 50:
        name = name2 = ''.join(random.choices(string.ascii_uppercase + string.digits, k=8))
        vars()[name] = pt.Genome(name="plasmid", length=450,
                            transcript_degradation_rate=1e-2,
                            transcript_degradation_rate_ext=1e-2,
                            rnase_speed=20,
                            rnase_footprint=10)
        eps = np.random.normal(mu, sigma)
        f_new = f_old * (1.0 + eps)
        new_pol_strength = random.randint(1e10, 5e10)
        if f_new > f_old:
            #Accepting mutation
            prom = vars()[name].add_promoter(name="p1", start=1, stop=10,
                                 interactions={"rnapol": new_pol_strength})
        else:
            #Calculate fitness of new mutation
            probability = calc_fitness(f_new, f_old, Ne)
            f_old = probability
            if probability > random.random():
                #Accepting mutation
                prom = vars()[name].add_promoter(name="p1", start=1, stop=10,
                                     interactions={"rnapol": new_pol_strength})

        #recreated_genome(prom, name2)
        #Recreated three genome sequence in order to simulate in loop
        sim = pt.Model(cell_volume=8e-16)
        sim.seed(34)
        sim.add_polymerase(name="rnapol", copy_number=4, speed=40, footprint=10)
        sim.add_ribosome(copy_number=100, speed=30, footprint=10)

        vars()[name] = pt.Genome(name="plasmid", length=450,
                            transcript_degradation_rate=1e-2,
                            transcript_degradation_rate_ext=1e-2,
                            rnase_speed=20,
                            rnase_footprint=10)
        prom
        vars()[name].add_terminator(name="t1", start=449, stop=450,
                               efficiency={"rnapol": 1.0})

        vars()[name].add_gene(name="proteinX", start=26, stop=148,
                         rbs_start=(26 - 15), rbs_stop=26, rbs_strength=1e7)

        vars()[name].add_gene(name="proteinY", start=26 + 150, stop=148 + 150,
                         rbs_start=(26 - 15 + 150), rbs_stop=26 + 150, rbs_strength=1e7)

        vars()[name].add_gene(name="proteinZ", start=165 + 150, stop=298 + 150,
                         rbs_start=(165 - 15 + 150), rbs_stop=165 + 150, rbs_strength=1e7)
        sim.register_genome(vars()[name])
        sim.simulate(time_limit=240, time_step=1,
                     output = "three_genes_rnase2.tsv")
        #nf = pandas.read_csv("three_genes_rnase2.tsv", sep="/t", engine="python")
        #check_files(df, nf)
        i+=1

#Calculates the fitness of the new mutation
def calc_fitness(variant_fit, orig_fit, Ne):
    xi = orig_fit
    xj = variant_fit
    N = Ne

    if xj == xi:
        return(1.0)
    if xj > xi:
        return(1.0)
    else:
        variant_fit = (xj / xi) ** (2 * N - 2)
        return variant_fit

def check_files(target_file, new_file):
    for line1 in new_file:
        for line2 in target_file:
            if line1 == line2:
                break

'''def sum_of_squares(target_file, new_file):
    for row in new_file.iterrows():
        for row2 in target_file.iterrows():
            if (row2[3] - row[3]) > 0:
                new_pol_strength += 1
            elif(row2[3] - row[3]) < 0:
                new_pol_strength -= 1
            elif(row2[3] - row[3]) == 0:
                check_files(target_file, new_file)
                break'''

'''def recreated_genome(promoter, name):
    sim = pt.Model(cell_volume=8e-16)
    sim.seed(34)
    sim.add_polymerase(name="rnapol", copy_number=4, speed=40, footprint=10)
    sim.add_ribosome(copy_number=100, speed=30, footprint=10)

    name = pt.Genome(name="plasmid", length=450,
                        transcript_degradation_rate=1e-2,
                        transcript_degradation_rate_ext=1e-2,
                        rnase_speed=20,
                        rnase_footprint=10)
    promoter
    name.add_terminator(name="t1", start=449, stop=450,
                           efficiency={"rnapol": 1.0})

    name.add_gene(name="proteinX", start=26, stop=148,
                     rbs_start=(26 - 15), rbs_stop=26, rbs_strength=1e7)

    name.add_gene(name="proteinY", start=26 + 150, stop=148 + 150,
                     rbs_start=(26 - 15 + 150), rbs_stop=26 + 150, rbs_strength=1e7)

    name.add_gene(name="proteinZ", start=165 + 150, stop=298 + 150,
                     rbs_start=(165 - 15 + 150), rbs_stop=165 + 150, rbs_strength=1e7)
    sim.register_genome(name)
    sim.simulate(time_limit=240, time_step=1,
                 output = "three_genes_rnase2.tsv")'''




if __name__ == '__main__':
    main()
