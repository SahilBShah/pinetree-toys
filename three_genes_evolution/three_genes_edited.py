import pinetree as pt

#Creating starting three genes sequence
sim = pt.Model(cell_volume=8e-16)
sim.seed(34)
sim.add_polymerase(name="rnapol", copy_number=4, speed=40, footprint=10)
sim.add_ribosome(copy_number=100, speed=30, footprint=10)

plasmid = pt.Genome(name="plasmid", length=450,
                    transcript_degradation_rate=1e-2,
                    transcript_degradation_rate_ext=1e-2,
                    rnase_speed=20,
                    rnase_footprint=10)

# plasmid = pt.Genome(name="plasmid", length=300)

plasmid.add_promoter(name="p1", start=1, stop=10,
                     interactions={"rnapol": 2e10})
plasmid.add_terminator(name="t1", start=449, stop=450,
                       efficiency={"rnapol": 1.0})

plasmid.add_gene(name="proteinX", start=26, stop=121,
                 rbs_start=11, rbs_stop=26, rbs_strength=1e7)

plasmid.add_gene(name="proteinY", start=159, stop=280,
                 rbs_start=144, rbs_stop=159, rbs_strength=1e7)
plasmid.add_rnase_site(start=133, stop=143)

plasmid.add_gene(name="proteinZ", start=319, stop=449,
                 rbs_start=303, rbs_stop=318, rbs_strength=1e7)
sim.register_genome(plasmid)
sim.simulate(time_limit=240, time_step=1,
             output = "three_genes_test_file.tsv")
