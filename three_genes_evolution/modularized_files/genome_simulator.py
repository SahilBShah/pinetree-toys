import pinetree as pt

def simulate_genome(genome_tracker, file_name):
    """
    Calls the PineTree simulator to run genome simulations and outputs a file containing all relevant information to be compared to target input file.
    """

    genome_tracker['start'] = genome_tracker['start'].astype(int)
    genome_tracker['stop'] = genome_tracker['stop'].astype(int)
    sim = pt.Model(cell_volume=8e-16)
    sim.seed(34)
    sim.add_polymerase(name="rnapol", copy_number=4, speed=40, footprint=10)
    sim.add_ribosome(copy_number=100, speed=30, footprint=10)

    plasmid = pt.Genome(name="plasmid", length=450,
                        transcript_degradation_rate=1e-2,
                        transcript_degradation_rate_ext=1e-2,
                        rnase_speed=20,
                        rnase_footprint=10)
    plasmid.add_promoter(name="p1", start=1, stop=10,
                         interactions={"rnapol": genome_tracker.loc[('promoter1'), ('current_strength')]})
    if genome_tracker.loc[('promoter2'), ('start')] > 0:
        plasmid.add_promoter(name="p2", start=genome_tracker.loc[('promoter2'), ('start')], stop=genome_tracker.loc[('promoter2'), ('stop')],
                             interactions={"rnapol": genome_tracker.loc[('promoter2'), ('current_strength')]})
    if genome_tracker.loc[('promoter3'), ('start')] > 0:
        plasmid.add_promoter(name="p3", start=genome_tracker.loc[('promoter3'), ('start')], stop=genome_tracker.loc[('promoter3'), ('stop')],
                             interactions={"rnapol": genome_tracker.loc[('promoter3'), ('current_strength')]})
    if genome_tracker.loc[('rnase1'), ('start')] > 0:
        plasmid.add_rnase_site(start=genome_tracker.loc[('rnase1'), ('start')], stop=genome_tracker.loc[('rnase1'), ('stop')])
    if genome_tracker.loc[('rnase2'), ('start')] > 0:
        plasmid.add_rnase_site(start=genome_tracker.loc[('rnase2'), ('start')], stop=genome_tracker.loc[('rnase2'), ('stop')])
    if genome_tracker.loc[('rnase3'), ('start')] > 0:
        plasmid.add_rnase_site(start=genome_tracker.loc[('rnase3'), ('start')], stop=genome_tracker.loc[('rnase3'), ('stop')])
    if genome_tracker.loc[('terminator1'), ('start')] > 0:
        plasmid.add_terminator(name="t1", start=genome_tracker.loc[('terminator1'), ('start')], stop=genome_tracker.loc[('terminator1'), ('stop')],
                               efficiency={"rnapol": genome_tracker.loc[('terminator1'), ('current_strength')]})
    if genome_tracker.loc[('terminator2'), ('start')] > 0:
        plasmid.add_terminator(name="t2", start=genome_tracker.loc[('terminator2'), ('start')], stop=genome_tracker.loc[('terminator2'), ('stop')],
                               efficiency={"rnapol": genome_tracker.loc[('terminator2'), ('current_strength')]})
    if genome_tracker.loc[('terminator3'), ('start')] > 0:
        plasmid.add_terminator(name="t3", start=genome_tracker.loc[('terminator3'), ('start')], stop=genome_tracker.loc[('terminator3'), ('stop')],
                               efficiency={"rnapol": genome_tracker.loc[('terminator3'), ('current_strength')]})
    plasmid.add_gene(name="proteinX", start=genome_tracker.loc[('geneX'), ('start')], stop=genome_tracker.loc[('geneX'), ('stop')],
                     rbs_start=(genome_tracker.loc[('geneX'), ('start')]-15), rbs_stop=genome_tracker.loc[('geneX'), ('start')], rbs_strength=1e7)
    plasmid.add_gene(name="proteinY", start=genome_tracker.loc[('geneY'), ('start')], stop=genome_tracker.loc[('geneY'), ('stop')],
                     rbs_start=(genome_tracker.loc[('geneY'), ('start')]-15), rbs_stop=genome_tracker.loc[('geneY'), ('start')], rbs_strength=1e7)
    plasmid.add_gene(name="proteinZ", start=genome_tracker.loc[('geneZ'), ('start')], stop=genome_tracker.loc[('geneZ'), ('stop')],
                     rbs_start=(genome_tracker.loc[('geneZ'), ('start')]-15), rbs_stop=genome_tracker.loc[('geneZ'), ('start')], rbs_strength=1e7)
    sim.register_genome(plasmid)
    sim.simulate(time_limit=240, time_step=1, output=file_name)
