import pinetree as pt
import pandas
import random
import numpy as np

#Taking in target tsv file
#df = pandas.read_csv(input("Enter .tsv file name: "), sep="/t", engine="python")
df = pandas.read_csv("three_genes_rnase.tsv", sep="/t", engine="python")

def main():
    f_old = 1.0
    sigma = 1
    eps = random.randint(0,5) / 100 * sigma
    Ne = 5
    i = 0
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

    prom = plasmid.add_promoter(name="p1", start=1, stop=10,
                         interactions={"rnapol": 4e10})
    plasmid.add_terminator(name="t1", start=449, stop=450,
                           efficiency={"rnapol": 1.0})

    plasmid.add_gene(name="proteinX", start=26, stop=148,
                     rbs_start=(26 - 15), rbs_stop=26, rbs_strength=1e7)

    plasmid.add_gene(name="proteinY", start=26 + 150, stop=148 + 150,
                     rbs_start=(26 - 15 + 150), rbs_stop=26 + 150, rbs_strength=1e7)

    plasmid.add_gene(name="proteinZ", start=165 + 150, stop=298 + 150,
                     rbs_start=(165 - 15 + 150), rbs_stop=165 + 150, rbs_strength=1e7)

    while i < 10:
        eps = random.randint(0,5) / 100 * sigma
        f_new = f_old * (1.0 + eps)
        new_pol_strength = random.randint(1e10, 5e10)
        #calculate fitness of new mutation
        probability = calc_fitness(f_new, f_old, Ne)
        f_old = probability
        print("Polymerase Strength = ", new_pol_strength)
        if f_new > f_old:
            prom = plasmid.add_promoter(name="p1", start=1, stop=10,
                                 interactions={"rnapol": new_pol_strength})
        else:
            probability = calc_fitness(f_new, f_old, Ne)
            f_old = probability
            if probability > random.random():
                prom = plasmid.add_promoter(name="p1", start=1, stop=10,
                                     interactions={"rnapol": new_pol_strength})
        #Add check equivalence bw tsv files feature
        print("f_new = ", f_new)
        print("f_old = ", f_old)
        prom
        sim.register_genome(plasmid)
        sim.simulate(time_limit=240, time_step=1,
                     output="three_genes_rnase2.tsv")
        nf = pandas.read_csv("three_genes_rnase2.tsv", sep="/t", engine="python")
        check_files(df, nf)
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




if __name__ == '__main__':
    main()
