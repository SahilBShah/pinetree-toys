import file_setup
import pandas as pd
# import random
# import yaml
# import save
# import numpy as np
# import csv

def genome_setup():
    """
    Sets up input files with necessary data for evolution to occur
    """

    input_df = input("Enter tsv file name: ")
    df = pd.read_csv(input_df, header=0, sep='\t')
    df = file_setup.rearrange_file(df)

    return df
