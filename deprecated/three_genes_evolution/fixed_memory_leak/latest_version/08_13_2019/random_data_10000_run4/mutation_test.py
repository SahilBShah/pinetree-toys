import sum_of_squares
import mutation_accepted
import save
import fitness_score
import random
from os_genome import call_pt
import yaml

def test_mutation(df, genome_tracker, i):
    """
    Mutation is tested to either be accepted or rejected based on calculated fitness values or based on a probability to allow for some neutral mutations to be accepted.
    """

    if i in [1, 5000, 10000]:
        genome_tracker['output_file_name'] = "gen_{}_data.tsv".format(i)
        with open('new_gene.yml', 'w') as f:
            yaml.dump(genome_tracker, f)
        genome_output_file_name = "genome_{}_tracker.yml".format(i)
        call_pt.pt_call()
        with open(genome_output_file_name, 'w') as f:
            yaml.dump(genome_tracker, f)
        f.close()

    if genome_tracker['f_new'] > genome_tracker['f_old']:
        mutation_accepted.accept_mutation(df, genome_tracker, i)
    else:
        #Calculate fitness of new mutation
        probability = fitness_score.calc_fitness(genome_tracker['f_new'], genome_tracker['f_old'])
        genome_tracker['f_old'] = probability
        with open('new_gene.yml', 'w') as f:
            yaml.dump(genome_tracker, f)
        f.close()
        if probability > random.random():
            mutation_accepted.accept_mutation(df, genome_tracker, i)
