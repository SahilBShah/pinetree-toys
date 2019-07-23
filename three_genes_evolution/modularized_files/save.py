import pandas as pd
import csv

def save_file(file, file_name):

    output_file = pd.DataFrame(file)
    output_file.to_csv(file_name, sep='\t')
