import pinetree as pt
import pandas as pd
import string
import random
import math
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
    i = 1
    gen = 0

    input_df = input("Enter tsv file name: ")
    df = pd.read_csv(input_df, header=0, sep='\t')
    df = rearrange_file(df)
    genome_tracker = pd.read_csv("gene_tracker.tsv", header=0, sep='\t')
    genome_tracker.reset_index()
    genome_tracker.set_index('species', inplace=True)
    genome_tracker.loc[('f_old'), ('value')] = 1.0

    #Setting random strengths for promoters and terminators
    genome_tracker.loc[('promoter1'), ('previous_strength')] = genome_tracker.loc[('promoter1'), ('current_strength')]
    genome_tracker.loc[('promoter1'), ('current_strength')] = random.randint(0, 3e10)

    possibilities = [modify_promoter(genome_tracker), modify_rnase(genome_tracker), modify_terminator(genome_tracker)]

    #Start of evolution program
    while i < 10001:

        random.choice(possibilities)

        eps = np.random.normal(mu, sigma)
        genome_tracker.loc[('f_new'), ('value')] = genome_tracker.loc[('f_old'), ('value')] * (1.0 + eps)
        test_mutation(df, genome_tracker, i)

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
    print(sos)
    return sos

#Removes unnecessary rows and columns in produced file
def rearrange_file(file):

    file = file[file['species'].isin(['proteinX', 'proteinY', 'proteinZ'])]
    file = file[['time', 'species', 'transcript']]
    file['time'] = file['time'].round().astype(int)
    file = file.pivot(index='time', columns='species', values='transcript')
    file = file.fillna(0.0)
    return file

def test_mutation(df, genome_tracker, i):

    f_new = genome_tracker.loc[('f_new'), ('value')]
    f_old = genome_tracker.loc[('f_old'), ('value')]

    if i in [1, 5000, 10000]:
        output_file_name = "gen_{}_data.tsv".format(i)
        genome_output_file_name = "genome_tracker_{}.tsv".format(i)
        genome_start.simulate_genome(genome_tracker, output_file_name)
        save_file(genome_tracker, genome_output_file_name)

    #f_new > f_old
    if f_new > f_old:
        accept_mutation(df, genome_tracker, i)
    else:
        #Calculate fitness of new mutation
        probability = calc_fitness(f_new, f_old, Ne)
        #f_old
        genome_tracker.loc[('f_old'), ('value')] = probability
        if probability > random.random():
            accept_mutation(df, genome_tracker, i)

def accept_mutation(df, genome_tracker, i):

    global gen
    global sos_list
    global all_sos_list
    output_file_name = "three_genes_replicated.tsv"
    genome_output_file_name = "genome_tracker_" + str(i) +".tsv"

    #Accepting mutation
    genome_start.simulate_genome(genome_tracker, output_file_name)
    #Taking in new file and removing unnecessary rows and columns
    nf = pd.read_csv("three_genes_replicated.tsv", header=0, sep='\t')
    nf = rearrange_file(nf)
    sos = calc_sum_of_squares(df, nf)
    all_sos_list.append(sos)
    #Accepts mutation only if sum of squares value decreases
    if sos <= sos_list[-1]:
        output_file_name = "three_genes_replicated_" + str(i) + ".tsv"
        save_file(nf, output_file_name)
        save_file(genome_tracker, genome_output_file_name)
        sos_list.append(sos)
        print(sos_list)
        gen+=1

def modify_promoter(genome_tracker):

    promoter_modification = ['add', 'remove']
    chosen_modification = random.choice(promoter_modification)
    promoter_slots = ['A', 'B']
    poly_sigma = 1e5

    if chosen_modification == 'add':
        promoter_possibilities = ["promoter1", "promoter2", "promoter3"]
        chosen_promoter = random.choice(promoter_possibilities)
        if chosen_promoter == 'promoter2':
            regionA = 'region2a'
            regionB = 'region2b'
            rnase = 'rnase1'
            terminator = 'terminator1'
            gene_endA = 123
            gene_endB = 134
        elif chosen_promoter == 'promoter3':
            regionA = 'region3a'
            regionB = 'region3b'
            rnase = 'rnase2'
            terminator = 'terminator2'
            gene_endA = 282
            gene_endB = 293

        if chosen_promoter != "promoter1":
            #Adding in a promoter between genes 1 and 2
            if (genome_tracker.loc[(rnase), ('start')] >= genome_tracker.loc[(regionA), ('start')]) and (genome_tracker.loc[(rnase), ('start')] <= genome_tracker.loc[(regionA), ('stop')]):
                promoter_slots.remove('A')
            if (genome_tracker.loc[(terminator), ('start')] >= genome_tracker.loc[(regionA), ('start')]) and (genome_tracker.loc[(terminator), ('start')] <= genome_tracker.loc[(regionA), ('stop')]):
                promoter_slots.remove('A')
            if (genome_tracker.loc[(rnase), ('start')] >= genome_tracker.loc[(regionB), ('start')]) and (genome_tracker.loc[(rnase), ('start')] <= genome_tracker.loc[(regionB), ('stop')]):
                promoter_slots.remove('B')
            if (genome_tracker.loc[(terminator), ('start')] >= genome_tracker.loc[(regionB), ('start')]) and (genome_tracker.loc[(terminator), ('stop')] <= genome_tracker.loc[(regionB), ('stop')]):
                promoter_slots.remove('B')
            if promoter_slots == []:
                genome_tracker.loc[(chosen_promoter), ('start')] = 0
                genome_tracker.loc[(chosen_promoter), ('stop')] = 0
            else:
                available_slot = random.choice(promoter_slots)
                if available_slot == 'A':
                    prom_start = genome_tracker.loc[(chosen_promoter), ('start')] = random.randint(genome_tracker.loc[(regionA), ('start')], gene_endA)
                    prom_stop = genome_tracker.loc[(chosen_promoter), ('stop')] = prom_start + 9
                if available_slot == 'B':
                    prom_start = genome_tracker.loc[(chosen_promoter), ('start')] = random.randint(genome_tracker.loc[(regionB), ('start')], gene_endB)
                    prom_stop = genome_tracker.loc[(chosen_promoter), ('stop')] = prom_start + 9

        #Determining polymerase strength
        poly_eps = np.random.normal(mu, poly_sigma)
        #Determines new polymerase strength to reduce sum of squares value
        pol_strength = genome_tracker.loc[(chosen_promoter), ('current_strength')] + poly_eps
        while pol_strength < 0:
            poly_eps = np.random.normal(mu, poly_sigma)
            pol_strength = genome_tracker.loc[(chosen_promoter), ('current_strength')] + poly_eps
            if pol_strength > 3e10:
                poly_eps = np.random.normal(mu, poly_sigma)
                pol_strength = genome_tracker.loc[(chosen_promoter), ('current_strength')] + poly_eps
        genome_tracker.loc[(chosen_promoter), ('previous_strength')] = genome_tracker.loc[(chosen_promoter), ('current_strength')]
        genome_tracker.loc[(chosen_promoter), ('current_strength')] = pol_strength

    if chosen_modification == 'remove':
        promoter_possibilities = ['promoter2', 'promoter3']
        chosen_promoter = random.choice(promoter_possibilities)
        #Removing promoters
        genome_tracker.loc[(chosen_promoter), ('start')] = 0
        genome_tracker.loc[(chosen_promoter), ('stop')] = 0

def modify_rnase(genome_tracker):

    rnase_modification = ['remove', 'add']
    chosen_modification = random.choice(rnase_modification)
    rnase_possibilities = ["rnase1", "rnase2", "rnase3"]
    chosen_rnase = random.choice(rnase_possibilities)
    rnase_slots = ['A', 'B', 'C']

    if chosen_modification == 'add':
        if chosen_rnase == 'rnase2':
            regionA = 'region2a'
            regionB = 'region2b'
            regionC = 'region2c'
            promoter = 'promoter2'
            terminator = 'terminator1'
            gene_endA = 122
            gene_endB = 133
            gene_endC = 149
        elif chosen_rnase == 'rnase3':
            regionA = 'region3a'
            regionB = 'region3b'
            regionC = 'region3c'
            promoter = 'promoter3'
            terminator = 'terminator2'
            gene_endA = 281
            gene_endB = 292
            gene_endC = 308

        if chosen_rnase == "rnase1":
            #Adds rnase after first promoter
            rnase1_start = genome_tracker.loc[('rnase1'), ('start')] = random.randint(genome_tracker.loc[('region1'), ('start')], 15)
            rnase1_stop = genome_tracker.loc[('rnase1'), ('stop')] = rnase1_start + 10
        else:
            #Adds rnase after second promoter
            if (genome_tracker.loc[(promoter), ('start')] >= genome_tracker.loc[(regionA), ('start')]) and (genome_tracker.loc[(promoter), ('start')] <= genome_tracker.loc[(regionA), ('stop')]):
                rnase_slots.remove('A')
            if (genome_tracker.loc[(terminator), ('start')] >= genome_tracker.loc[(regionA), ('start')]) and (genome_tracker.loc[(terminator), ('start')] <= genome_tracker.loc[(regionA), ('stop')]):
                rnase_slots.remove('A')
            if (genome_tracker.loc[(promoter), ('start')] >= genome_tracker.loc[(regionB), ('start')]) and (genome_tracker.loc[(promoter), ('start')] <= genome_tracker.loc[(regionB), ('stop')]):
                rnase_slots.remove('B')
            if (genome_tracker.loc[(terminator), ('start')] >= genome_tracker.loc[(regionB), ('start')]) and (genome_tracker.loc[(terminator), ('start')] <= genome_tracker.loc[(regionB), ('stop')]):
                rnase_slots.remove('B')
            if (genome_tracker.loc[(terminator), ('start')] >= genome_tracker.loc[(regionC), ('start')]) and (genome_tracker.loc[(terminator), ('start')] <= genome_tracker.loc[(regionC), ('stop')]):
                rnase_slots.remove('C')
            if rnase_slots == []:
                genome_tracker.loc[(chosen_rnase), ('start')] = 0
                genome_tracker.loc[(chosen_rnase), ('stop')] = 0
            else:
                available_slot = random.choice(rnase_slots)
                if available_slot == 'A':
                    rnase_start = genome_tracker.loc[(chosen_rnase), ('start')] = random.randint(genome_tracker.loc[(regionA), ('start')], gene_endA)
                    rnase_stop = genome_tracker.loc[(chosen_rnase), ('stop')] = rnase_start + 10
                if available_slot == 'B':
                    rnase_start = genome_tracker.loc[(chosen_rnase), ('start')] = random.randint(genome_tracker.loc[(regionB), ('start')], gene_endB)
                    rnase_stop = genome_tracker.loc[(chosen_rnase), ('stop')] = rnase_start + 10
                if available_slot == 'C':
                    rnase_start = genome_tracker.loc[(chosen_rnase), ('start')] = random.randint(genome_tracker.loc[(regionC), ('start')], gene_endC)
                    rnase_stop = genome_tracker.loc[(chosen_rnase), ('stop')] = rnase_start + 10

    if chosen_modification == 'remove':
        genome_tracker.loc[(chosen_rnase), ('start')] = 0
        genome_tracker.loc[(chosen_rnase), ('stop')] = 0

def modify_terminator(genome_tracker):

    terminator_modification = ['remove', 'add']
    chosen_modification = random.choice(terminator_modification)
    terminator_possibilities = ["terminator1", "terminator2", "terminator3"]
    chosen_terminator = random.choice(terminator_possibilities)
    terminator_slots = ['A', 'B', 'C']
    term_sigma = 1.0

    if chosen_modification == 'add':
        if chosen_terminator == 'terminator1':
            regionA = 'region2a'
            regionB = 'region2b'
            regionC = 'region2c'
            rnase = 'rnase1'
            promoter = 'promoter2'
            gene_endA = 131
            gene_endB = 142
            gene_endC = 158
        elif chosen_terminator == 'terminator2':
            regionA = 'region3a'
            regionB = 'region3b'
            regionC = 'region3c'
            rnase = 'rnase2'
            promoter = 'promoter3'
            gene_endA = 290
            gene_endB = 301
            gene_endC = 317

        if chosen_terminator == "terminator3":
            #Adds terminator after third gene
            genome_tracker.loc[(chosen_terminator), ('start')] = 449
            genome_tracker.loc[(chosen_terminator), ('stop')] = 450
        else:
            #Adds terminator after first gene
            if (genome_tracker.loc[(promoter), ('start')] >= genome_tracker.loc[(regionA), ('start')]) and (genome_tracker.loc[(promoter), ('start')] <= genome_tracker.loc[(regionA), ('stop')]):
                terminator_slots.remove('A')
            if (genome_tracker.loc[(rnase), ('start')] >= genome_tracker.loc[(regionA), ('start')]) and (genome_tracker.loc[(rnase), ('start')] <= genome_tracker.loc[(regionA), ('stop')]):
                terminator_slots.remove('A')
            if (genome_tracker.loc[(promoter), ('start')] >= genome_tracker.loc[(regionB), ('start')]) and (genome_tracker.loc[(promoter), ('start')] <= genome_tracker.loc[(regionB), ('stop')]):
                terminator_slots.remove('B')
            if (genome_tracker.loc[(rnase), ('start')] >= genome_tracker.loc[(regionB), ('start')]) and (genome_tracker.loc[(rnase), ('start')] <= genome_tracker.loc[(regionB), ('stop')]):
                terminator_slots.remove('B')
            if (genome_tracker.loc[(rnase), ('start')] >= genome_tracker.loc[(regionC), ('start')]) and (genome_tracker.loc[(rnase), ('start')] <= genome_tracker.loc[(regionC), ('stop')]):
                terminator_slots.remove('C')
            if terminator_slots == []:
                genome_tracker.loc[(chosen_terminator), ('start')] = 0
                genome_tracker.loc[(chosen_terminator), ('stop')] = 0
            else:
                available_slot = random.choice(terminator_slots)
                if available_slot == 'A':
                    term_start = genome_tracker.loc[(chosen_terminator), ('start')] = random.randint(genome_tracker.loc[(regionA), ('start')], gene_endA)
                    term_stop = genome_tracker.loc[(chosen_terminator), ('stop')] = term_start + 1
                if available_slot == 'B':
                    term_start = genome_tracker.loc[(chosen_terminator), ('start')] = random.randint(genome_tracker.loc[(regionB), ('start')], gene_endB)
                    term_stop = genome_tracker.loc[(chosen_terminator), ('stop')] = term_start + 1
                if available_slot == 'C':
                    term_start = genome_tracker.loc[(chosen_terminator), ('start')] = random.randint(genome_tracker.loc[(regionC), ('start')], gene_endC)
                    term_stop = genome_tracker.loc[(chosen_terminator), ('stop')] = term_start + 1

        #Determining terminator polymerase efficiency rate
        term_eps = np.random.normal(mu, term_sigma)
        #Determines new terminator efficiency value to reduce sum of squares value
        term_efficiency = genome_tracker.loc[(chosen_terminator), ('current_strength')] + term_eps
        while term_efficiency < 0 or term_efficiency > 1:
            term_eps = np.random.normal(mu, term_sigma)
            term_efficiency = genome_tracker.loc[(chosen_terminator), ('current_strength')] + term_eps
        genome_tracker.loc[(chosen_terminator), ('current_strength')] = term_efficiency
        genome_tracker.loc[(chosen_terminator), ('previous_strength')] = genome_tracker.loc[(chosen_terminator), ('current_strength')]

    if terminator_modification == 'remove':
        genome_tracker.loc[(chosen_terminator), ('start')] = 0
        genome_tracker.loc[(chosen_terminator), ('stop')] = 0
        genome_tracker.loc[(chosen_terminator), ('previous_strength')] = genome_tracker.loc[(chosen_terminator), ('current_strength')]
        genome_tracker.loc[(chosen_terminator), ('current_strength')] = 0.0

def save_file(file, file_name):

    output_file = pd.DataFrame(file)
    output_file.to_csv(file_name, sep='\t')

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
                         interactions={"rnapol": genome_tracker.loc[('promoter1'), ('current_strength')]})
    if genome_tracker.loc[('promoter2'), ('start')] > 0:
        plasmid.add_promoter(name="p2", start=genome_tracker.loc[('promoter2'), ('start')], stop=genome_tracker.loc[('promoter2'), ('stop')],
                             interactions={"rnapol": genome_tracker.loc[('promoter2'), ('current_strength')]})
    if genome_tracker.loc[('promoter3'), ('start')] > 0:
        plasmid.add_promoter(name="p3", start=genome_tracker.loc[('promoter3'), ('start')], stop=genome_tracker.loc[('promoter3'), ('stop')],
                             interactions={"rnapol": genome_tracker.loc[('promoter3'), ('current_strength')]})
    if genome_tracker.loc[('rnase1'), ('start')] > 0:
        plasmid.add_rnase_site(start=genome_tracker.loc[('rnase1'), ('start')], stop=genome_tracker.loc[('rnase1'), ('stop')])
    if genome_tracker.loc[('rnase2'), ('start')] > 0:
        plasmid.add_rnase_site(start=genome_tracker.loc[('rnase2'), ('start')], stop=genome_tracker.loc[('rnase2'), ('stop')])
    if genome_tracker.loc[('rnase3'), ('start')] > 0:
        plasmid.add_rnase_site(start=genome_tracker.loc[('rnase3'), ('start')], stop=genome_tracker.loc[('rnase3'), ('stop')])
    if genome_tracker.loc[('terminator1'), ('start')] > 0:
        plasmid.add_terminator(name="t1", start=genome_tracker.loc[('terminator1'), ('start')], stop=genome_tracker.loc[('terminator1'), ('stop')],
                               efficiency={"rnapol": genome_tracker.loc[('terminator1'), ('current_strength')]})
    if genome_tracker.loc[('terminator2'), ('start')] > 0:
        plasmid.add_terminator(name="t2", start=genome_tracker.loc[('terminator2'), ('start')], stop=genome_tracker.loc[('terminator2'), ('stop')],
                               efficiency={"rnapol": genome_tracker.loc[('terminator2'), ('current_strength')]})
    if genome_tracker.loc[('terminator3'), ('start')] > 0:
        plasmid.add_terminator(name="t3", start=genome_tracker.loc[('terminator3'), ('start')], stop=genome_tracker.loc[('terminator3'), ('stop')],
                               efficiency={"rnapol": genome_tracker.loc[('terminator3'), ('current_strength')]})
    plasmid.add_gene(name="proteinX", start=genome_tracker.loc[('geneX'), ('start')], stop=genome_tracker.loc[('geneX'), ('stop')],
                     rbs_start=(genome_tracker.loc[('geneX'), ('start')]-15), rbs_stop=genome_tracker.loc[('geneX'), ('start')], rbs_strength=1e7)
    plasmid.add_gene(name="proteinY", start=genome_tracker.loc[('geneY'), ('start')], stop=genome_tracker.loc[('geneY'), ('stop')],
                     rbs_start=(genome_tracker.loc[('geneY'), ('start')]-15), rbs_stop=genome_tracker.loc[('geneY'), ('start')], rbs_strength=1e7)
    plasmid.add_gene(name="proteinZ", start=genome_tracker.loc[('geneZ'), ('start')], stop=genome_tracker.loc[('geneZ'), ('stop')],
                     rbs_start=(genome_tracker.loc[('geneZ'), ('start')]-15), rbs_stop=genome_tracker.loc[('geneZ'), ('start')], rbs_strength=1e7)
    sim.register_genome(plasmid)
    sim.simulate(time_limit=240, time_step=1, output=file_name)

if __name__ == '__main__':
    main()
