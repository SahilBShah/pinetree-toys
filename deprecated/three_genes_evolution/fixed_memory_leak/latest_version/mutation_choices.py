import numpy as np
import random
import pandas as pd
import yaml

global mu
mu = 0.0

def modify_promoter(genome_tracker_new):
    """
    Promoters are either added with randomized polymerase strengths or removed from the genome all together.
    """

    promoter_modification = ['add', 'remove', 'modify']
    chosen_prom_modification = random.choice(promoter_modification)
    promoter_slots = ['A', 'B']
    poly_sigma = 1e5

    if chosen_prom_modification == 'add':
        promoter_possibilities = ["promoter1", "promoter2", "promoter3"]
        chosen_promoter = random.choice(promoter_possibilities)
        if chosen_promoter == 'promoter2':
            regionA = 'region2a'
            regionB = 'region2b'
            rnase = 'rnase1'
            terminator = 'terminator1'
            gene_endA = 122
            gene_endB = 133
        elif chosen_promoter == 'promoter3':
            regionA = 'region3a'
            regionB = 'region3b'
            rnase = 'rnase2'
            terminator = 'terminator2'
            gene_endA = 281
            gene_endB = 292

        if chosen_promoter != "promoter1":
            #Adding in a promoter between genes 1 and 2
            if (genome_tracker_new[rnase]['start'] >= genome_tracker_new[regionA]['start']) and (genome_tracker_new[rnase]['start'] <= genome_tracker_new[regionA]['stop']):
                promoter_slots.remove('A')
            if (genome_tracker_new[terminator]['start'] >= genome_tracker_new[regionA]['start']) and (genome_tracker_new[terminator]['start'] <= genome_tracker_new[regionA]['stop']):
                promoter_slots.remove('A')
            if (genome_tracker_new[rnase]['start'] >= genome_tracker_new[regionB]['start']) and (genome_tracker_new[rnase]['start'] <= genome_tracker_new[regionB]['stop']):
                promoter_slots.remove('B')
            if (genome_tracker_new[terminator]['start'] >= genome_tracker_new[regionB]['start']) and (genome_tracker_new[terminator]['stop'] <= genome_tracker_new[regionB]['stop']):
                promoter_slots.remove('B')
            if promoter_slots == []:
                genome_tracker_new[chosen_promoter]['start'] = 0
                genome_tracker_new[chosen_promoter]['stop'] = 0
            else:
                available_slot = random.choice(promoter_slots)
                if available_slot == 'A':
                    prom_start = genome_tracker_new[chosen_promoter]['start'] = random.randint(genome_tracker_new[regionA]['start'], gene_endA)
                    prom_stop = genome_tracker_new[chosen_promoter]['stop'] = prom_start + 9
                if available_slot == 'B':
                    prom_start = genome_tracker_new[chosen_promoter]['start'] = random.randint(genome_tracker_new[regionB]['start'], gene_endB)
                    prom_stop = genome_tracker_new[chosen_promoter]['stop'] = prom_start + 9

        genome_tracker_new[chosen_promoter]['previous_strength'] = genome_tracker_new[chosen_promoter]['current_strength']
        genome_tracker_new[chosen_promoter]['current_strength'] = 10e4

    if chosen_prom_modification == "modify":
        print("promoter modify")
        promoter_possibilities = ['promoter1']
        if genome_tracker_new['promoter2']['start'] > 0:
            promoter_possibilities.append('promoter2')
        if genome_tracker_new['promoter3']['start'] > 0:
            promoter_possibilities.append('promoter3')
        chosen_promoter = random.choice(promoter_possibilities)
        genome_tracker_new[chosen_promoter]['previous_strength'] = genome_tracker_new[chosen_promoter]['current_strength']
        prom_strength = genome_tracker_new[chosen_promoter]['current_strength'] * np.random.normal(1, 0.1)
        while 3e10 < prom_strength < 0.0:
            prom_strength = genome_tracker_new[chosen_promoter]['current_strength'] * np.random.normal(1, 0.1)
        genome_tracker_new[chosen_promoter]['current_strength'] = prom_strength

    if chosen_prom_modification == 'remove':
        promoter_possibilities = ['promoter2', 'promoter3']
        chosen_promoter = random.choice(promoter_possibilities)
        #Removing promoters
        genome_tracker_new[chosen_promoter]['start'] = 0
        genome_tracker_new[chosen_promoter]['stop'] = 0

    with open('new_gene.yml', 'w') as f:
        yaml.dump(genome_tracker_new, f, default_flow_style=False)

def modify_rnase(genome_tracker_new):
    """
    Rnases are added or removed.
    """

    rnase_modification = ['remove', 'add']
    chosen_rnase_modification = random.choice(rnase_modification)
    rnase_possibilities = ["rnase1", "rnase2", "rnase3"]
    chosen_rnase = random.choice(rnase_possibilities)
    rnase_slots = ['A', 'B', 'C']

    if chosen_rnase_modification == 'add':
        if chosen_rnase == 'rnase2':
            regionA = 'region2a'
            regionB = 'region2b'
            regionC = 'region2c'
            promoter = 'promoter2'
            terminator = 'terminator1'
            gene_endA = 122
            gene_endB = 133
            gene_endC = 144
        elif chosen_rnase == 'rnase3':
            regionA = 'region3a'
            regionB = 'region3b'
            regionC = 'region3c'
            promoter = 'promoter3'
            terminator = 'terminator2'
            gene_endA = 281
            gene_endB = 292
            gene_endC = 303

        if chosen_rnase == "rnase1":
            #Adds rnase after first promoter
            rnase1_start = genome_tracker_new['rnase1']['start'] = random.randint(genome_tracker_new['region1']['start'], 15)
            rnase1_stop = genome_tracker_new['rnase1']['stop'] = rnase1_start + 10
        else:
            #Adds rnase after second promoter
            if (genome_tracker_new[promoter]['start'] >= genome_tracker_new[regionA]['start']) and (genome_tracker_new[promoter]['start'] <= genome_tracker_new[regionA]['stop']):
                rnase_slots.remove('A')
            if (genome_tracker_new[terminator]['start'] >= genome_tracker_new[regionA]['start']) and (genome_tracker_new[terminator]['start'] <= genome_tracker_new[regionA]['stop']):
                rnase_slots.remove('A')
            if (genome_tracker_new[promoter]['start'] >= genome_tracker_new[regionB]['start']) and (genome_tracker_new[promoter]['start'] <= genome_tracker_new[regionB]['stop']):
                rnase_slots.remove('B')
            if (genome_tracker_new[terminator]['start'] >= genome_tracker_new[regionB]['start']) and (genome_tracker_new[terminator]['start'] <= genome_tracker_new[regionB]['stop']):
                rnase_slots.remove('B')
            if (genome_tracker_new[terminator]['start'] >= genome_tracker_new[regionC]['start']) and (genome_tracker_new[terminator]['start'] <= genome_tracker_new[regionC]['stop']):
                rnase_slots.remove('C')
            if rnase_slots == []:
                genome_tracker_new[chosen_rnase]['start'] = 0
                genome_tracker_new[chosen_rnase]['stop'] = 0
            else:
                available_slot = random.choice(rnase_slots)
                if available_slot == 'A':
                    rnase_start = genome_tracker_new[chosen_rnase]['start'] = random.randint(genome_tracker_new[regionA]['start'], gene_endA)
                    rnase_stop = genome_tracker_new[chosen_rnase]['stop'] = rnase_start + 10
                if available_slot == 'B':
                    rnase_start = genome_tracker_new[chosen_rnase]['start'] = random.randint(genome_tracker_new[regionB]['start'], gene_endB)
                    rnase_stop = genome_tracker_new[chosen_rnase]['stop'] = rnase_start + 10
                if available_slot == 'C':
                    rnase_start = genome_tracker_new[chosen_rnase]['start'] = random.randint(genome_tracker_new[regionC]['start'], gene_endC)
                    rnase_stop = genome_tracker_new[chosen_rnase]['stop'] = rnase_start + 10

    if chosen_rnase_modification == 'remove':
        genome_tracker_new[chosen_rnase]['start'] = 0
        genome_tracker_new[chosen_rnase]['stop'] = 0

    with open('new_gene.yml', 'w') as f:
        yaml.dump(genome_tracker_new, f, default_flow_style=False)

def modify_terminator(genome_tracker_new):
    """
    Terminators are either added with randomized terminator efficiencies or removed all together.
    """

    terminator_modification = ['remove', 'add', 'modify']
    chosen_term_modification = random.choice(terminator_modification)
    terminator_possibilities = ["terminator1", "terminator2", "terminator3"]
    chosen_terminator = random.choice(terminator_possibilities)
    terminator_slots = ['A', 'B', 'C']
    term_sigma = 1.0

    if chosen_term_modification == 'add':
        if chosen_terminator == 'terminator1':
            regionA = 'region2a'
            regionB = 'region2b'
            regionC = 'region2c'
            rnase = 'rnase1'
            promoter = 'promoter2'
            gene_endA = 131
            gene_endB = 142
            gene_endC = 158
        elif chosen_terminator == 'terminator2':
            regionA = 'region3a'
            regionB = 'region3b'
            regionC = 'region3c'
            rnase = 'rnase2'
            promoter = 'promoter3'
            gene_endA = 290
            gene_endB = 301
            gene_endC = 317

        if chosen_terminator == "terminator3":
            #Adds terminator after third gene
            genome_tracker_new[chosen_terminator]['start'] = 449
            genome_tracker_new[chosen_terminator]['stop'] = 450
        else:
            #Adds terminator after first gene
            if (genome_tracker_new[promoter]['start'] >= genome_tracker_new[regionA]['start']) and (genome_tracker_new[promoter]['start'] <= genome_tracker_new[regionA]['stop']):
                terminator_slots.remove('A')
            if (genome_tracker_new[rnase]['start'] >= genome_tracker_new[regionA]['start']) and (genome_tracker_new[rnase]['start'] <= genome_tracker_new[regionA]['stop']):
                terminator_slots.remove('A')
            if (genome_tracker_new[promoter]['start'] >= genome_tracker_new[regionB]['start']) and (genome_tracker_new[promoter]['start'] <= genome_tracker_new[regionB]['stop']):
                terminator_slots.remove('B')
            if (genome_tracker_new[rnase]['start'] >= genome_tracker_new[regionB]['start']) and (genome_tracker_new[rnase]['start'] <= genome_tracker_new[regionB]['stop']):
                terminator_slots.remove('B')
            if (genome_tracker_new[rnase]['start'] >= genome_tracker_new[regionC]['start']) and (genome_tracker_new[rnase]['start'] <= genome_tracker_new[regionC]['stop']):
                terminator_slots.remove('C')
            if terminator_slots == []:
                genome_tracker_new[chosen_terminator]['start'] = 0
                genome_tracker_new[chosen_terminator]['stop'] = 0
            else:
                available_slot = random.choice(terminator_slots)
                if available_slot == 'A':
                    term_start = genome_tracker_new[chosen_terminator]['start'] = random.randint(genome_tracker_new[regionA]['start'], gene_endA)
                    term_stop = genome_tracker_new[chosen_terminator]['stop'] = term_start + 1
                if available_slot == 'B':
                    term_start = genome_tracker_new[chosen_terminator]['start'] = random.randint(genome_tracker_new[regionB]['start'], gene_endB)
                    term_stop = genome_tracker_new[chosen_terminator]['stop'] = term_start + 1
                if available_slot == 'C':
                    term_start = genome_tracker_new[chosen_terminator]['start'] = random.randint(genome_tracker_new[regionC]['start'], gene_endC)
                    term_stop = genome_tracker_new[chosen_terminator]['stop'] = term_start + 1

        genome_tracker_new[chosen_terminator]['previous_strength'] = genome_tracker_new[chosen_terminator]['current_strength']
        genome_tracker_new[chosen_terminator]['current_strength'] = 0.85

    if chosen_term_modification == "modify":
        print("terminator modify")
        terminator_possibilities = []
        if genome_tracker_new['terminator1']['start'] > 0:
            terminator_possibilities.append('terminator1')
        if genome_tracker_new['terminator2']['start'] > 0:
            terminator_possibilities.append('terminator2')
        if genome_tracker_new['terminator3']['start'] > 0:
            terminator_possibilities.append('terminator3')
        if terminator_possibilities != []:
            chosen_terminator = random.choice(terminator_possibilities)
            genome_tracker_new[chosen_terminator]['previous_strength'] = genome_tracker_new[chosen_terminator]['current_strength']
            term_eff = genome_tracker_new[chosen_terminator]['current_strength'] * np.random.normal(1, 0.1)
            while 1.0 < term_eff < 0.0:
                term_eff = genome_tracker_new[chosen_terminator]['current_strength'] * np.random.normal(1, 0.1)
            genome_tracker_new[chosen_terminator]['current_strength'] = term_eff

    if chosen_term_modification == 'remove':
        genome_tracker_new[chosen_terminator]['start'] = 0
        genome_tracker_new[chosen_terminator]['stop'] = 0
        genome_tracker_new[chosen_terminator]['previous_strength'] = genome_tracker_new[chosen_terminator]['current_strength']
        genome_tracker_new[chosen_terminator]['current_strength'] = 0.0

    with open('new_gene.yml', 'w') as f:
        yaml.dump(genome_tracker_new, f, default_flow_style=False)
