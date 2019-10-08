import sum_of_squares
from os_genome import call_pt
import pandas as pd
import file_setup

def test_mutation(df):
    """
    Mutation is tested to either be accepted or rejected based on calculated fitness values or based on a probability to allow for some neutral mutations to be accepted.
    """

    call_pt.pt_call()
    saveable_nf = nf = pd.read_csv("three_genes_replicated.tsv", header=0, sep='\t')
    nf = file_setup.rearrange_file(nf)
    return(sum_of_squares.calc_sum_of_squares(df, nf))
