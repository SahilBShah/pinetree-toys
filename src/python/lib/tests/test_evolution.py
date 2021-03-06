#Common imports
import copy
import datetime
import filecmp
from glob import glob
import pandas as pd
import os
import pytest
import yaml

os.chdir(os.path.dirname(os.path.abspath(__file__)))

import sys

sys.path.insert(1, '../../')

#Evolution package
from ...evolution import *


def test_enumerate_mutation_options():

    #Verify the mutation options are correctly assigned
    with open(os.path.dirname(os.path.abspath(__file__))+'/inputs/testing.yml', 'r') as gene_parameters:
        genome_tracker = yaml.safe_load(gene_parameters)
    mutation_possibilities = enumerate_mutation_options(genome_tracker)
    assert mutation_possibilities != {}
    assert 'promoter_2.remove' in mutation_possibilities
    assert mutation_possibilities['promoter_2.remove'] == 'remove'
    assert mutation_possibilities['promoter_2.remove'] != 'add'
    assert 'terminator_1.add' in mutation_possibilities
    assert mutation_possibilities['terminator_1.add'] == 'add'
    assert 'rnase_2.modify0' in mutation_possibilities
    assert mutation_possibilities['rnase_2.modify0'] == 'modify'

def test_evolution():

    #Change directory to src
    os.chdir('../../')
    #Run evolutionary program
    os.system('python evolution.py paper_data1_arrange1.tsv starting.yml 10000 10 3')
    #Check to see if evolutionary run completed
    assert os.path.isfile('../../results/paper_data1_arrange1/rep10000/final/rmse_data.tsv')
    assert os.path.isfile('../../results/paper_data1_arrange1/rep10000/final/expression_pattern_best.tsv')
    assert os.path.isfile('../../results/paper_data1_arrange1/rep10000/final/gene_best.yml')
    #Check if starting genome yaml file was created
    assert os.path.isfile('../../results/paper_data1_arrange1/rep10000/config.yml')
    assert os.path.isfile('../../results/paper_data1_arrange1/rep10000/gene_0.yml')
    assert os.path.isfile('../../results/paper_data1_arrange1/rep10000/expression_pattern_0.tsv')
    #Check if expression data and genome files are saved for each accepted mutation
    assert len(glob('../../results/paper_data1_arrange1/rep10000/expression_pattern_*.tsv')) != 0
    assert len(glob('../../results/paper_data1_arrange1/rep10000/gene_*.yml'))
    os.system('rm -rf ../../results/paper_data1_arrange1')
    os.chdir(os.getcwd())
