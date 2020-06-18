#Common imports
import copy
import numpy as np
import pandas as pd
import random
from scipy import stats
import yaml

#lib imports
import mutation_analysis

#promoter_starting_strength = 10e6
# max_promoter_strength = 10e13
# min_promoter_strength = 10e5
#promoter_offset = 34
#promoter_min_space = 35
#rnase_starting_strength = 5e-3
# max_rnase_strength = 1.0
# min_rnase_strength = 0.0
#rnase_offset = 9
#rnase_min_space = 10
# max_terminator_strength = 1.0
# min_terminator_strength = 0.0
#terminator_starting_strength = 0.3
#terminator_offset = 29
#terminator_min_space = 30

#{element: [genome_shift, starting_strengths]}
element_dict = {'promoter': [34, 10e6], 'rnase': [9, 5e-3], 'terminator': [29, 0.3]}
#{element: [min_strength, max_strength]}
element_strengths_range = {'promoter': [10e5, 10e13], 'rnase': [0.0, 1.0], 'terminator': [0.0, 1.0]}


def add_element(genome_tracker_new, output_dir, num_genes, deg_rate, element_choice):
    """
    Proposes a muttion that adds an element on to the genome.
    Input(s):
    genome_tracker_new is the dataframe containing the most recent edited genomic data.
    output_dir is the path to the directory in which all the saved files are stored by the program.
    num_genes refers to the number of genes in the genome.
    deg_rate is a command line argument that specifies if rnase degredation rates should be individually specified or not.
    element_choice is the element selected to be added.
    """

    genome_elements = []
    index_list = []
    spaces_dict = {}

    region_choice = 'region_{}'.format(element_choice.split('_')[1])
    #Appends elements that are already on the genome on to a list
    promoter = 'promoter_{}'.format(region_choice.split('_')[1])
    terminator = 'terminator_{}'.format(region_choice.split('_')[1])
    rnase = 'rnase_{}'.format(region_choice.split('_')[1])
    gene = 'gene_{}'.format(region_choice.split('_')[1])

    #Appends base pair values in which an element can be placed between
    if int(region_choice.split('_')[1]) != num_genes:
        #If a promoter is present on the genome
        if genome_tracker_new[promoter]['start'] > 0:
            #If the region selected is not the region before the first gene
            if region_choice != 'region_0':
                genome_elements.append(genome_tracker_new[promoter]['start'])
            genome_elements.append(genome_tracker_new[promoter]['stop'])
        #If an rnase is present on the genome
        if genome_tracker_new[rnase]['start'] > 0:
            genome_elements.append(genome_tracker_new[rnase]['start'])
            genome_elements.append(genome_tracker_new[rnase]['stop'])
    #If the region selected is not the region before the first gene
    if region_choice != 'region_0':
        genome_elements.append(genome_tracker_new[region_choice]['start'])
        #If a terminator is present on the genome
        if genome_tracker_new[terminator]['start'] > 0:
            genome_elements.append(genome_tracker_new[terminator]['start'])
            #If the end of the terminator is not the end of the genome
            if genome_tracker_new[terminator]['stop'] != genome_tracker_new['length_of_genome']:
                genome_elements.append(genome_tracker_new[terminator]['stop'])
    #If the last terminator ends at the end of the genome then only append the end of the genome to avoid duplicate values
    if int(region_choice.split('_')[1]) == num_genes and genome_tracker_new[terminator]['stop'] != genome_tracker_new['length_of_genome']:
        genome_elements.append(genome_tracker_new['length_of_genome'])
    genome_elements.append(genome_tracker_new[region_choice]['stop'])

    #Determines areas within the intergenic regions that are available
    space_index = 0
    highest_position = max(genome_elements)
    while True:
        #Add positions between elements to a dictionary
        starting_position = min(genome_elements)
        genome_elements.remove(starting_position)
        ending_position = min(genome_elements)
        genome_elements.remove(ending_position)
        spaces_dict['space{}'.format(space_index)] = dict(start=starting_position, stop=ending_position)
        if ending_position == highest_position or genome_elements == []:
            break
        space_index+=1
    #Identifies the area in which to add an element in
    key = random.choice(list(spaces_dict.keys()))
    starting_position = spaces_dict[key]['start']
    ending_position = spaces_dict[key]['stop']

    #Base pair value to shift the genome by depending on the element chosen
    genome_shift = element_dict[element_choice.split('_')[0]][0] + 1
    if starting_position + 1 == ending_position:
        ending_position+=1
    #Adds the element to the genome with its starting strength
    start = genome_tracker_new[element_choice]['start'] = random.randint(starting_position+1, ending_position)
    stop = genome_tracker_new[element_choice]['stop'] = start + element_dict[element_choice.split('_')[0]][0]
    genome_tracker_new[element_choice]['previous_strength'] = genome_tracker_new[element_choice]['current_strength']
    genome_tracker_new[element_choice]['current_strength'] = element_dict[element_choice.split('_')[0]][1]
    #Increase the length of the genome by the genome_shift value
    genome_tracker_new = expand_genome(genome_tracker_new, num_genes, int(region_choice.split('_')[1]), genome_shift, element_choice)

    return genome_tracker_new

def remove_element(genome_tracker_new, output_dir, num_genes, deg_rate, element_choice):
    """
    Proposes a mutation that removes an element from the genome.
    Input(s):
    genome_tracker_new is the dataframe containing the most recent edited genomic data.
    output_dir is the path to the directory in which all the saved files are stored by the program.
    num_genes refers to the number of genes in the genome.
    deg_rate is a command line argument that specifies if rnase degredation rates should be individually specified or not.
    element_choice is the element selected to be removed.
    """

    #Decrease the length of the genome by the genome_shift value
    genome_tracker_new = shrink_genome(genome_tracker_new, num_genes, int(element_choice.split('_')[1]), element_dict[element_choice.split('_')[0]][0]+1, element_choice)
    #Removes the selected element from the genome
    genome_tracker_new[element_choice]['start'] = 0
    genome_tracker_new[element_choice]['stop'] = 0
    genome_tracker_new[element_choice]['previous_strength'] = genome_tracker_new[element_choice]['current_strength']
    genome_tracker_new[element_choice]['current_strength'] = 0

    return genome_tracker_new

def modify_element(genome_tracker_new, output_dir, num_genes, deg_rate, element_choice):
    """
    Modify an element's (that is present on the genome) strength.
    Input(s):
    genome_tracker_new is the dataframe containing the most recent edited genomic data.
    output_dir is the path to the directory in which all the saved files are stored by the program.
    num_genes refers to the number of genes in the genome.
    deg_rate is a command line argument that specifies if rnase degredation rates should be individually specified or not.
    element_choice is the element selected to be modified.
    """

    #Modify the strength of the selected element present on the genome
    genome_tracker_new[element_choice]['previous_strength'] = genome_tracker_new[element_choice]['current_strength']
    genome_tracker_new[element_choice]['current_strength'] = genome_tracker_new[element_choice]['current_strength'] * np.random.normal(1, 0.2)
    #If the modified strength goes below the minimum threshold or above the maximum threshold, continue to modify the element's strength
    while genome_tracker_new[element_choice]['current_strength'] < element_strengths_range[element_choice.split('_')[0]][0] or genome_tracker_new[element_choice]['current_strength'] > element_strengths_range[element_choice.split('_')[0]][1]:
        genome_tracker_new[element_choice]['current_strength'] = genome_tracker_new[element_choice]['previous_strength'] * np.random.normal(1, 0.2)

    return genome_tracker_new

def expand_genome(genome_tracker_new, num_genes, region_choice, genome_shift, element_choice):
    """
    When an element is added, the genome length increases.
    Input(s):
    genome_tracker_new is the dataframe containing the most recent edited genomic data.
    num_genes refers to the number of genes in the genome.
    region_choice is the region the selected element resides in.
    genome_shift is the amount of base pairs to add to the genome when an element is added in.
    element_choice is the element selected to be added.
    """

    #Increases the genome size if an element is added
    for beg_point in range(region_choice, num_genes+1):

        #List possible elements on the genome
        region = 'region_{}'.format(beg_point)
        terminator = 'terminator_{}'.format(beg_point)
        promoter = 'promoter_{}'.format(beg_point)
        rnase = 'rnase_{}'.format(beg_point)
        gene_offsetted = 'gene_{}'.format(beg_point+1)

        #Increase the end of the selected region
        genome_tracker_new[region]['stop']+=genome_shift
        #If the iteration is not the current region or region before the first gene, increase the beginning position of the curent region
        if beg_point != region_choice and region != 'region_0':
            genome_tracker_new[region]['start']+=genome_shift
        #If the iteration is not the region before the first gene or the region after the last gene
        if region != 'region_0' and beg_point != num_genes:
            #If the selected element is not a promoter
            if element_choice != promoter:
                #If the promoter come after the selected element and is present on the genome, increase the position of the promoter
                if genome_tracker_new[promoter]['start'] >= genome_tracker_new[element_choice]['start'] and genome_tracker_new[promoter]['start'] > 0:
                    genome_tracker_new[promoter]['start']+=genome_shift
                    genome_tracker_new[promoter]['stop']+=genome_shift
        #If the iteration is not the region before the first gene
        if region != 'region_0':
            #If the selected element is not a terminator
            if element_choice != terminator:
                #If the terminator come after the selected element and is present on the genome, increase the position of the terminator
                if genome_tracker_new[terminator]['start'] >= genome_tracker_new[element_choice]['start'] and genome_tracker_new[terminator]['start'] > 0:
                    genome_tracker_new[terminator]['start']+=genome_shift
                    genome_tracker_new[terminator]['stop']+=genome_shift
        #If the iteration is not the region after the last gene, increase the position of all the genes after the selected element
        if beg_point != num_genes:
            genome_tracker_new[gene_offsetted]['start']+=genome_shift
            genome_tracker_new[gene_offsetted]['stop']+= genome_shift
            #If the selected element is not an rnase
            if element_choice != rnase:
                #If the rnase come after the selected element and is present on the genome, increase the position of the rnase
                if genome_tracker_new[rnase]['start'] >= genome_tracker_new[element_choice]['start'] and genome_tracker_new[rnase]['start'] > 0:
                    genome_tracker_new[rnase]['start']+=genome_shift
                    genome_tracker_new[rnase]['stop']+=genome_shift
    #Increase the overall length of the genome
    genome_tracker_new['length_of_genome']+=genome_shift

    return genome_tracker_new

def shrink_genome(genome_tracker_new, num_genes, region_choice, genome_shift, element_choice):
    """
    When an element is removed, the genome length decreases.
    Input(s):
    genome_tracker_new is the dataframe containing the most recent edited genomic data.
    num_genes refers to the number of genes in the genome.
    region_choice is the region the selected element resides in.
    genome_shift is the amount of base pairs to remove from the genome when an element is removed.
    element_choice is the element selected to be removed.
    """

    #Decreases the genome size if an element is added
    for beg_point in range(region_choice, num_genes+1):

        #List possible elements on the genome
        region = 'region_{}'.format(beg_point)
        terminator = 'terminator_{}'.format(beg_point)
        promoter = 'promoter_{}'.format(beg_point)
        rnase = 'rnase_{}'.format(beg_point)
        gene_offsetted = 'gene_{}'.format(beg_point+1)

        #Decrease the end of the selected region
        genome_tracker_new[region]['stop']-=genome_shift
        #If the iteration is not the current region or region before the first gene, decrease the beginning position of the curent region
        if beg_point != region_choice and region != 'region_0':
            genome_tracker_new[region]['start']-=genome_shift
        #If the iteration is not the region before the first gene or the region after the last gene
        if region != 'region_0' and beg_point != num_genes:
            #If the selected element is not a promoter
            if element_choice != promoter:
                #If the promoter come after the selected element and is present on the genome, decrease the position of the promoter
                if genome_tracker_new[promoter]['start'] >= genome_tracker_new[element_choice]['start'] and genome_tracker_new[promoter]['start'] > 0:
                    genome_tracker_new[promoter]['start']-=genome_shift
                    genome_tracker_new[promoter]['stop']-=genome_shift
        #If the iteration is not the region before the first gene
        if region != 'region_0':
            #If the selected element is not a terminator
            if element_choice != terminator:
                #If the terminator come after the selected element and is present on the genome, decrease the position of the terminator
                if genome_tracker_new[terminator]['start'] >= genome_tracker_new[element_choice]['start'] and genome_tracker_new[terminator]['start'] > 0:
                    genome_tracker_new[terminator]['start']-=genome_shift
                    genome_tracker_new[terminator]['stop']-=genome_shift
        #If the iteration is not the region after the last gene, decrease the position of all the genes after the selected element
        if beg_point != num_genes:
            genome_tracker_new[gene_offsetted]['start']-=genome_shift
            genome_tracker_new[gene_offsetted]['stop']-= genome_shift
            #If the selected element is not an rnase
            if element_choice != rnase:
                #If the rnase come after the selected element and is present on the genome, decrease the position of the rnase
                if genome_tracker_new[rnase]['start'] >= genome_tracker_new[element_choice]['start'] and genome_tracker_new[rnase]['start'] > 0:
                    genome_tracker_new[rnase]['start']-=genome_shift
                    genome_tracker_new[rnase]['stop']-=genome_shift
    #Decrease the overall length of the genome
    genome_tracker_new['length_of_genome']-=genome_shift

    return genome_tracker_new

def cleanup_genome(output_dir, target_file, sse_df, num_genes, deg_rate):
    """
    Removes elements that have low strengths that are statistically insignificant to the overal gene expression pattern produced from the best found architecture.
    Input(s):
    output_dir is the path to the directory in which all the saved files are stored by the program.
    target_file is the user-inputted tsv file containing transcript abundances for each gene.
    sse_df is the dataframe containing all the sum of squared error data at each generation.
    num_genes refers to the number of genes in the genome.
    deg_rate is a command line argument that specifies if rnase degredation rates should be individually specified or not.
    """

    remove_elements = []

    #Get index with lowest sum of squared error value
    sse_df = sse_df[sse_df['Accepted'] == 'yes']
    min_sse_df = sse_df[sse_df.SSE == sse_df.SSE.min()]
    min_sse_index = min_sse_df.iloc[-1]['Iteration']

    #Open genome architecture file
    with open(output_dir+'gene_{}.yml'.format(min_sse_index), 'r') as gene_elements:
        genome_tracker_best = yaml.safe_load(gene_elements)
    genome_tracker_saved = copy.deepcopy(genome_tracker_best)
    #Calculate SSE for best genome and get a range of values
    ss_best = mutation_analysis.analyze_mutation(genome_tracker_best, output_dir, target_file, 20, deg_rate, True)

    #Iterate through each element and remove them and compare SSE values of new architecture to the best architecture previously found
    for gene in range(num_genes+1):
        promoter = 'promoter_{}'.format(gene)
        terminator = 'terminator_{}'.format(gene)
        rnase = 'rnase_{}'.format(gene)
        if gene != 0:
            if genome_tracker_best[terminator]['start'] > 0:
                genome_tracker_saved = remove_element(genome_tracker_saved, output_dir, num_genes, deg_rate, terminator)
                #Get SSE values range to compare to the best found architecture
                ss_comp = mutation_analysis.analyze_mutation(genome_tracker_saved, output_dir, target_file, 20, deg_rate, True)
                if stats.ttest_ind(ss_best, ss_comp)[1] >= 0.05:
                    remove_elements.append((terminator, stats.ttest_ind(ss_best, ss_comp)[1]))
                #Set the genome architecture file back to original state
                genome_tracker_saved = copy.deepcopy(genome_tracker_best)
        if gene != num_genes and gene != 0:
            if genome_tracker_best[promoter]['start'] > 0:
                genome_tracker_saved = remove_element(genome_tracker_saved, output_dir, num_genes, deg_rate, promoter)
                #Get SSE values range to compare to the best found architecture
                ss_comp = mutation_analysis.analyze_mutation(genome_tracker_saved, output_dir, target_file, 20, deg_rate, True)
                if stats.ttest_ind(ss_best, ss_comp)[1] >= 0.05:
                    remove_elements.append((promoter, stats.ttest_ind(ss_best, ss_comp)[1]))
                #Set the genome architecture file back to original state
                genome_tracker_saved = copy.deepcopy(genome_tracker_best)
        if gene != num_genes:
            if genome_tracker_best[rnase]['start'] > 0:
                genome_tracker_saved = remove_element(genome_tracker_saved, output_dir, num_genes, deg_rate, rnase)
                #Get SSE values range to compare to the best found architecture
                ss_comp = mutation_analysis.analyze_mutation(genome_tracker_saved, output_dir, target_file, 20, deg_rate, True)
                if stats.ttest_ind(ss_best, ss_comp)[1] >= 0.05:
                    remove_elements.append((rnase, stats.ttest_ind(ss_best, ss_comp)[1]))
                #Set the genome architecture file back to original state
                genome_tracker_saved = copy.deepcopy(genome_tracker_best)



    #Remove all elements that did not significantly alter the gene expression pattern produced
    remove_elements.sort(key=lambda x: (-x[1],x[0]))
    for element in remove_elements:
        genome_tracker_saved = remove_element(genome_tracker_saved, output_dir, num_genes, deg_rate, element[0])
        ss_comp = mutation_analysis.analyze_mutation(genome_tracker_saved, output_dir, target_file, 20, deg_rate, True)
        if stats.ttest_ind(ss_best, ss_comp)[1] >= 0.05:
            genome_tracker_best = copy.deepcopy(genome_tracker_saved)
        else:
            break

    #Saves best genome architecture found
    with open(output_dir+'final/gene_best.yml', 'w') as save_yaml:
        yaml.dump(genome_tracker_best, save_yaml)
    mutation_analysis.analyze_mutation(genome_tracker_saved, output_dir, target_file, 20, deg_rate)
    save_df = pd.read_csv(output_dir+"expression_pattern.tsv", header=0, sep='\t')
    save_df.to_csv(output_dir+"final/expression_pattern_best.tsv", sep='\t', index=False)

    return
