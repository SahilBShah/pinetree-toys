import sum_of_squares
import mutation_accepted
import save
import fitness_score
import random
from os_genome import call_pt
import yaml
import pandas as pd
import file_setup

def test_mutation(df, genome_tracker, i):
    """
    Mutation is tested to either be accepted or rejected based on calculated fitness values or based on a probability to allow for some neutral mutations to be accepted.
    """

    '''if i in [1, 5000, 10000]:
        genome_tracker['output_file_name'] = "gen_{}_data.tsv".format(i)
        with open('new_gene.yml', 'w') as f:
            yaml.dump(genome_tracker, f)
        genome_output_file_name = "genome_{}_tracker.yml".format(i)
        call_pt.pt_call()
        with open(genome_output_file_name, 'w') as f:
            yaml.dump(genome_tracker, f)
        f.close()'''

    genome_tracker['output_file_name'] = "three_genes_replicated.tsv"
    with open('new_gene.yml', 'w') as f:
        yaml.dump(genome_tracker, f)
    f.close()

    call_pt.pt_call()
    saveable_nf = nf = pd.read_csv("three_genes_replicated.tsv", header=0, sep='\t')
    nf = file_setup.rearrange_file(nf)
    #if i != 1:
    return(sum_of_squares.calc_sum_of_squares(df, nf))
    #Taking in new file and removing unnecessary rows and columns
    # saveable_nf = nf = pd.read_csv("three_genes_replicated.tsv", header=0, sep='\t')
    # nf = file_setup.rearrange_file(nf)
    # if i != 1:
    #     genome_tracker['f_new'] = sum_of_squares.calc_sum_of_squares(df, nf)
    #     print('f_new:', genome_tracker['f_new'])
    #     with open('new_gene.yml', 'w') as f:
    #         yaml.dump(genome_tracker, f)
    #     f.close()
    # #Calculate fitness of new mutation
    # probability = fitness_score.calc_fitness(genome_tracker['f_new'], genome_tracker['f_old'])
    # print('probability:', probability)
    # if probability > random.random():
    #     genome_tracker['f_old'] = genome_tracker['f_new']
    #     with open('new_gene.yml', 'w') as f:
    #         yaml.dump(genome_tracker, f)
    #     f.close()
    #     mutation_accepted.accept_mutation(df, genome_tracker, i, saveable_nf, nf)
