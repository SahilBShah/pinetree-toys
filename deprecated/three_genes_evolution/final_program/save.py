import pandas as pd
import csv

def save_file(file, file_name):
    """
    Saves files as needed.
    """

    file.to_csv(file_name, sep='\t', index=False)
