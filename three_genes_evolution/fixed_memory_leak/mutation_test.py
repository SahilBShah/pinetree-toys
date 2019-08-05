import sum_of_squares
import mutation_accepted
import genome_simulator
import save
import fitness_score
import random

def test_mutation(df, genome_tracker, i):
    """
    Mutation is tested to either be accepted or rejected based on calculated fitness values or based on a probability to allow for some neutral mutations to be accepted.
    """

    if i in [1, 5000, 10000]:
        output_file_name = "gen_{}_data.tsv".format(i)
        genome_output_file_name = "genome_{}_tracker.tsv".format(i)
        genome_simulator.simulate_genome(genome_tracker, output_file_name)
        save.save_file(genome_tracker, genome_output_file_name)

    if genome_tracker['f_new'] > genome_tracker['f_old']:
        mutation_accepted.accept_mutation(df, genome_tracker, i)
    else:
        #Calculate fitness of new mutation
        probability = fitness_score.calc_fitness(genome_tracker['f_new'], genome_tracker['f_old'])
        genome_tracker['f_old'] = probability
        if probability > random.random():
            mutation_accepted.accept_mutation(df, genome_tracker, i)
