import numpy as np
import random
import pandas as pd

global mu
mu = 0.0

def modify_promoter(genome_tracker):
    """
    Promoters are either added with randomized polymerase strengths or removed from the genome all together.
    """

    promoter_modification = ['add', 'remove']
    chosen_modification = random.choice(promoter_modification)
    promoter_slots = ['A', 'B']
    poly_sigma = 1e5

    if chosen_modification == 'add':
        promoter_possibilities = ["promoter1", "promoter2", "promoter3"]
        chosen_promoter = random.choice(promoter_possibilities)
        if chosen_promoter == 'promoter2':
            regionA = 'region2a'
            regionB = 'region2b'
            rnase = 'rnase1'
            terminator = 'terminator1'
            gene_endA = 123
            gene_endB = 134
        elif chosen_promoter == 'promoter3':
            regionA = 'region3a'
            regionB = 'region3b'
            rnase = 'rnase2'
            terminator = 'terminator2'
            gene_endA = 282
            gene_endB = 293

        if chosen_promoter != "promoter1":
            #Adding in a promoter between genes 1 and 2
            if (genome_tracker.loc[(rnase), ('start')] >= genome_tracker.loc[(regionA), ('start')]) and (genome_tracker.loc[(rnase), ('start')] <= genome_tracker.loc[(regionA), ('stop')]):
                promoter_slots.remove('A')
            if (genome_tracker.loc[(terminator), ('start')] >= genome_tracker.loc[(regionA), ('start')]) and (genome_tracker.loc[(terminator), ('start')] <= genome_tracker.loc[(regionA), ('stop')]):
                promoter_slots.remove('A')
            if (genome_tracker.loc[(rnase), ('start')] >= genome_tracker.loc[(regionB), ('start')]) and (genome_tracker.loc[(rnase), ('start')] <= genome_tracker.loc[(regionB), ('stop')]):
                promoter_slots.remove('B')
            if (genome_tracker.loc[(terminator), ('start')] >= genome_tracker.loc[(regionB), ('start')]) and (genome_tracker.loc[(terminator), ('stop')] <= genome_tracker.loc[(regionB), ('stop')]):
                promoter_slots.remove('B')
            if promoter_slots == []:
                genome_tracker.loc[(chosen_promoter), ('start')] = 0
                genome_tracker.loc[(chosen_promoter), ('stop')] = 0
            else:
                available_slot = random.choice(promoter_slots)
                if available_slot == 'A':
                    prom_start = genome_tracker.loc[(chosen_promoter), ('start')] = random.randint(genome_tracker.loc[(regionA), ('start')], gene_endA)
                    prom_stop = genome_tracker.loc[(chosen_promoter), ('stop')] = prom_start + 9
                if available_slot == 'B':
                    prom_start = genome_tracker.loc[(chosen_promoter), ('start')] = random.randint(genome_tracker.loc[(regionB), ('start')], gene_endB)
                    prom_stop = genome_tracker.loc[(chosen_promoter), ('stop')] = prom_start + 9

        #Determining polymerase strength
        poly_eps = np.random.normal(mu, poly_sigma)
        #Determines new polymerase strength to reduce sum of squares value
        pol_strength = genome_tracker.loc[(chosen_promoter), ('current_strength')] + poly_eps
        while pol_strength < 0:
            poly_eps = np.random.normal(mu, poly_sigma)
            pol_strength = genome_tracker.loc[(chosen_promoter), ('current_strength')] + poly_eps
            if pol_strength > 3e10:
                poly_eps = np.random.normal(mu, poly_sigma)
                pol_strength = genome_tracker.loc[(chosen_promoter), ('current_strength')] + poly_eps
        genome_tracker.loc[(chosen_promoter), ('previous_strength')] = genome_tracker.loc[(chosen_promoter), ('current_strength')]
        genome_tracker.loc[(chosen_promoter), ('current_strength')] = pol_strength

    if chosen_modification == 'remove':
        promoter_possibilities = ['promoter2', 'promoter3']
        chosen_promoter = random.choice(promoter_possibilities)
        #Removing promoters
        genome_tracker.loc[(chosen_promoter), ('start')] = 0
        genome_tracker.loc[(chosen_promoter), ('stop')] = 0

def modify_rnase(genome_tracker):
    """
    Rnases are added or removed.
    """

    rnase_modification = ['remove', 'add']
    chosen_modification = random.choice(rnase_modification)
    rnase_possibilities = ["rnase1", "rnase2", "rnase3"]
    chosen_rnase = random.choice(rnase_possibilities)
    rnase_slots = ['A', 'B', 'C']

    if chosen_modification == 'add':
        if chosen_rnase == 'rnase2':
            regionA = 'region2a'
            regionB = 'region2b'
            regionC = 'region2c'
            promoter = 'promoter2'
            terminator = 'terminator1'
            gene_endA = 122
            gene_endB = 133
            gene_endC = 149
        elif chosen_rnase == 'rnase3':
            regionA = 'region3a'
            regionB = 'region3b'
            regionC = 'region3c'
            promoter = 'promoter3'
            terminator = 'terminator2'
            gene_endA = 281
            gene_endB = 292
            gene_endC = 308

        if chosen_rnase == "rnase1":
            #Adds rnase after first promoter
            rnase1_start = genome_tracker.loc[('rnase1'), ('start')] = random.randint(genome_tracker.loc[('region1'), ('start')], 15)
            rnase1_stop = genome_tracker.loc[('rnase1'), ('stop')] = rnase1_start + 10
        else:
            #Adds rnase after second promoter
            if (genome_tracker.loc[(promoter), ('start')] >= genome_tracker.loc[(regionA), ('start')]) and (genome_tracker.loc[(promoter), ('start')] <= genome_tracker.loc[(regionA), ('stop')]):
                rnase_slots.remove('A')
            if (genome_tracker.loc[(terminator), ('start')] >= genome_tracker.loc[(regionA), ('start')]) and (genome_tracker.loc[(terminator), ('start')] <= genome_tracker.loc[(regionA), ('stop')]):
                rnase_slots.remove('A')
            if (genome_tracker.loc[(promoter), ('start')] >= genome_tracker.loc[(regionB), ('start')]) and (genome_tracker.loc[(promoter), ('start')] <= genome_tracker.loc[(regionB), ('stop')]):
                rnase_slots.remove('B')
            if (genome_tracker.loc[(terminator), ('start')] >= genome_tracker.loc[(regionB), ('start')]) and (genome_tracker.loc[(terminator), ('start')] <= genome_tracker.loc[(regionB), ('stop')]):
                rnase_slots.remove('B')
            if (genome_tracker.loc[(terminator), ('start')] >= genome_tracker.loc[(regionC), ('start')]) and (genome_tracker.loc[(terminator), ('start')] <= genome_tracker.loc[(regionC), ('stop')]):
                rnase_slots.remove('C')
            if rnase_slots == []:
                genome_tracker.loc[(chosen_rnase), ('start')] = 0
                genome_tracker.loc[(chosen_rnase), ('stop')] = 0
            else:
                available_slot = random.choice(rnase_slots)
                if available_slot == 'A':
                    rnase_start = genome_tracker.loc[(chosen_rnase), ('start')] = random.randint(genome_tracker.loc[(regionA), ('start')], gene_endA)
                    rnase_stop = genome_tracker.loc[(chosen_rnase), ('stop')] = rnase_start + 10
                if available_slot == 'B':
                    rnase_start = genome_tracker.loc[(chosen_rnase), ('start')] = random.randint(genome_tracker.loc[(regionB), ('start')], gene_endB)
                    rnase_stop = genome_tracker.loc[(chosen_rnase), ('stop')] = rnase_start + 10
                if available_slot == 'C':
                    rnase_start = genome_tracker.loc[(chosen_rnase), ('start')] = random.randint(genome_tracker.loc[(regionC), ('start')], gene_endC)
                    rnase_stop = genome_tracker.loc[(chosen_rnase), ('stop')] = rnase_start + 10

    if chosen_modification == 'remove':
        genome_tracker.loc[(chosen_rnase), ('start')] = 0
        genome_tracker.loc[(chosen_rnase), ('stop')] = 0

def modify_terminator(genome_tracker):
    """
    Terminators are either added with randomized terminator efficiencies or removed all together.
    """

    terminator_modification = ['remove', 'add']
    chosen_modification = random.choice(terminator_modification)
    terminator_possibilities = ["terminator1", "terminator2", "terminator3"]
    chosen_terminator = random.choice(terminator_possibilities)
    terminator_slots = ['A', 'B', 'C']
    term_sigma = 1.0

    if chosen_modification == 'add':
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
            genome_tracker.loc[(chosen_terminator), ('start')] = 449
            genome_tracker.loc[(chosen_terminator), ('stop')] = 450
        else:
            #Adds terminator after first gene
            if (genome_tracker.loc[(promoter), ('start')] >= genome_tracker.loc[(regionA), ('start')]) and (genome_tracker.loc[(promoter), ('start')] <= genome_tracker.loc[(regionA), ('stop')]):
                terminator_slots.remove('A')
            if (genome_tracker.loc[(rnase), ('start')] >= genome_tracker.loc[(regionA), ('start')]) and (genome_tracker.loc[(rnase), ('start')] <= genome_tracker.loc[(regionA), ('stop')]):
                terminator_slots.remove('A')
            if (genome_tracker.loc[(promoter), ('start')] >= genome_tracker.loc[(regionB), ('start')]) and (genome_tracker.loc[(promoter), ('start')] <= genome_tracker.loc[(regionB), ('stop')]):
                terminator_slots.remove('B')
            if (genome_tracker.loc[(rnase), ('start')] >= genome_tracker.loc[(regionB), ('start')]) and (genome_tracker.loc[(rnase), ('start')] <= genome_tracker.loc[(regionB), ('stop')]):
                terminator_slots.remove('B')
            if (genome_tracker.loc[(rnase), ('start')] >= genome_tracker.loc[(regionC), ('start')]) and (genome_tracker.loc[(rnase), ('start')] <= genome_tracker.loc[(regionC), ('stop')]):
                terminator_slots.remove('C')
            if terminator_slots == []:
                genome_tracker.loc[(chosen_terminator), ('start')] = 0
                genome_tracker.loc[(chosen_terminator), ('stop')] = 0
            else:
                available_slot = random.choice(terminator_slots)
                if available_slot == 'A':
                    term_start = genome_tracker.loc[(chosen_terminator), ('start')] = random.randint(genome_tracker.loc[(regionA), ('start')], gene_endA)
                    term_stop = genome_tracker.loc[(chosen_terminator), ('stop')] = term_start + 1
                if available_slot == 'B':
                    term_start = genome_tracker.loc[(chosen_terminator), ('start')] = random.randint(genome_tracker.loc[(regionB), ('start')], gene_endB)
                    term_stop = genome_tracker.loc[(chosen_terminator), ('stop')] = term_start + 1
                if available_slot == 'C':
                    term_start = genome_tracker.loc[(chosen_terminator), ('start')] = random.randint(genome_tracker.loc[(regionC), ('start')], gene_endC)
                    term_stop = genome_tracker.loc[(chosen_terminator), ('stop')] = term_start + 1

        #Determining terminator polymerase efficiency rate
        term_eps = np.random.normal(mu, term_sigma)
        #Determines new terminator efficiency value to reduce sum of squares value
        term_efficiency = genome_tracker.loc[(chosen_terminator), ('current_strength')] + term_eps
        while term_efficiency < 0 or term_efficiency > 1:
            term_eps = np.random.normal(mu, term_sigma)
            term_efficiency = genome_tracker.loc[(chosen_terminator), ('current_strength')] + term_eps
        genome_tracker.loc[(chosen_terminator), ('current_strength')] = term_efficiency
        genome_tracker.loc[(chosen_terminator), ('previous_strength')] = genome_tracker.loc[(chosen_terminator), ('current_strength')]

    if terminator_modification == 'remove':
        genome_tracker.loc[(chosen_terminator), ('start')] = 0
        genome_tracker.loc[(chosen_terminator), ('stop')] = 0
        genome_tracker.loc[(chosen_terminator), ('previous_strength')] = genome_tracker.loc[(chosen_terminator), ('current_strength')]
        genome_tracker.loc[(chosen_terminator), ('current_strength')] = 0.0
