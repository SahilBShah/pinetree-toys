import pandas as pd

#Removes unnecessary rows and columns in produced file
def rearrange_file(file):

    file = file[file['species'].isin(['proteinX', 'proteinY', 'proteinZ'])]
    file = file[['time', 'species', 'transcript']]
    file['time'] = file['time'].round().astype(int)
    file = file.pivot(index='time', columns='species', values='transcript')
    file = file.fillna(0.0)
    return file
