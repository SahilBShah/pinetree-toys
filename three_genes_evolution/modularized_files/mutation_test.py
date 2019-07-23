import sum_of_squares
import mutation_accepted
import genome_simulator
import save
import fitness_score
import random

def test_mutation(df, genome_tracker, i):

    f_new = genome_tracker.loc[('f_new'), ('value')]
    f_old = genome_tracker.loc[('f_old'), ('value')]

    if i in [1, 5000, 10000]:
        output_file_name = "gen_{}_data.tsv".format(i)
        genome_output_file_name = "genome_tracker_{}.tsv".format(i)
        genome_simulator.simulate_genome(genome_tracker, output_file_name)
        save.save_file(genome_tracker, genome_output_file_name)

    #f_new > f_old
    if f_new > f_old:
        mutation_accepted.accept_mutation(df, genome_tracker, i)
    else:
        #Calculate fitness of new mutation
        probability = fitness_score.calc_fitness(f_new, f_old)
        #f_old
        genome_tracker.loc[('f_old'), ('value')] = probability
        if probability > random.random():
            mutation_accepted.accept_mutation(df, genome_tracker, i)
