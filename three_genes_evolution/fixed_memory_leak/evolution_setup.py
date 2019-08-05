import pandas as pd
import csv
import file_setup
import random
import yaml

def genome_setup():
    """
    Sets up input files with necessary data for evolution to occur
    """

    input_df = input("Enter tsv file name: ")
    df = pd.read_csv(input_df, header=0, sep='\t')
    df = file_setup.rearrange_file(df)
    with open('gene.yml') as f:
        genome_tracker = yaml.safe_load(f)
    #Setting random strengths for promoters and terminators
    genome_tracker['promoter1']['previous_strength'] = genome_tracker['promoter1']['current_strength']
    genome_tracker['promoter1']['current_strength'] = random.randint(0, 3e10)
    #Setting arbitrary old fitness value
    genome_tracker['f_old'] = 1.0
    genome_tracker.close()

    return [df, genome_tracker]
