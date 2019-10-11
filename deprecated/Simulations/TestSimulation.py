import pinetree as pt

model = pt.Model(cell_volume=8e-16)
#Define genome's name & length
plasmid = pt.Genome(name="myplasmid", length=300)
#Add promoters, terminators, & genes (any order)
plasmid.add_promoter(name="phi1", start=1, stop=10, interactions={"rnapol": 2e8})
plasmid.add_terminator(name="t1", start=299, stop=300, efficiency={"rnapol": 1.0})
plasmid.add_gene(name="rnapol", start=26, stop=225, rbs_start=11, rbs_stop=26, rbs_strength=1e7)
plasmid.add_gene(name="proteinX", start=241, stop=280, rbs_start=226, rbs_stop=241, rbs_strength=1e7)
#Register genome with the model
model.register_genome(plasmid)
#Add polymerases THEN ribosomes
model.add_polymerase(name="rnapol", speed=40, footprint=10, copy_number=10)
model.add_ribosome(speed=30, footprint=10, copy_number=100)
#Species reactions bw molecular species
model.add_reaction(reactants=['proteinX', 'rnapol'],
                   products=['rnapol-X'],
                   rate_constant=1e-7)
#Run the simulation
model.simulate(time_limit=60, time_step=1, output="simulation.tsv")
