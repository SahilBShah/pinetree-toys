import pandas as pd

#Removes unnecessary rows and columns in produced file
def rearrange_file(file, genome_tracker_new):
    """
    Rearranges files so that the target input file can be compared to the simulated output file so that sum of squares can be easily calculated.
    Input(s):
    file is the transcript abundances file outputted from pinetree.
    genome_tracker_new is the dataframe containing the most recent edited genomic data.
    """

    protein_species = []

    #File is rearranged to have columns represent each gene's transcript abundances and each row represent the time
    for gene_number in range(1, genome_tracker_new['num_genes']+1):
        protein_species.append('protein{}'.format(gene_number))
    file = file[file['species'].isin(protein_species)]
    file = file[['time', 'species', 'transcript']]
    #Time is rounded so that files can be compared to one another
    file['time'] = file['time'].round().astype(int)
    file = file.pivot(index='time', columns='species', values='transcript')
    file = file.fillna(0.0)
    #If a protein is not present at a certain time period then no (0) transcripts are produced
    for protein in protein_species:
        if protein not in file.columns:
            file[protein] = 0.0
    file = file[protein_species]

    return file
