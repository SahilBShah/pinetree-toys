import pinetree as pt
import pandas

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
#plasmid.add_promoter(name="p1", start=122, stop=132,
                     #interactions={"rnapol": 2e10})
#plasmid.add_promoter(name="p1", start=292, stop=302,
                     #interactions={"rnapol": 2e10})
#plasmid.add_terminator(name="t1", start=449, stop=450,
                       #efficiency={"rnapol": 1.0})

plasmid.add_gene(name="proteinX", start=26, stop=121,
                 rbs_start=11, rbs_stop=26, rbs_strength=1e7)

plasmid.add_gene(name="proteinY", start=159, stop=280,
                 rbs_start=144, rbs_stop=159, rbs_strength=1e7)
#plasmid.add_rnase_site(start=133, stop=143)
#plasmid.add_rnase_site(start=281, stop=291)
plasmid.add_gene(name="proteinZ", start=319, stop=449,
                 rbs_start=303, rbs_stop=318, rbs_strength=1e7)
sim.register_genome(plasmid)
sim.simulate(time_limit=240, time_step=1,
             output = "three_genes_test_file.tsv")

'''gene1_list = [26, 121]
gene2_list = [159, 280]
gene3_list = [319, 449]
prom1_list = [122, 132]
prom2_list = [292, 302]
term1_list = [0, 0]
term2_list = [0, 0]
term3_list = [449, 450]
rnase1_list = [0, 0]
rnase2_list = [133, 143]
rnase3_list = [281, 291]

gene1_dataframe = pandas.DataFrame(gene1_list, columns=["Gene_Location"])
export_csv = gene1_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/target_genome/gene1_data.tsv", index=False)
gene2_dataframe = pandas.DataFrame(gene2_list, columns=["Gene_Location"])
export_csv = gene2_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/target_genome/gene2_data.tsv", index=False)
gene3_dataframe = pandas.DataFrame(gene3_list, columns=["Gene_Location"])
export_csv = gene3_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/target_genome/gene3_data.tsv", index=False)
prom1_dataframe = pandas.DataFrame(data=prom1_list, columns=["Promoter_Location"])
export_csv = prom1_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/target_genome/promoter1_data.tsv", index=False)
prom2_dataframe = pandas.DataFrame(data=prom2_list, columns=["Promoter_Location"])
export_csv = prom2_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/target_genome/promoter2_data.tsv", index=False)
term1_dataframe = pandas.DataFrame(term1_list, columns=["Terminator_Location"])
export_csv = term1_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/target_genome/term1_data.tsv", index=False)
term2_dataframe = pandas.DataFrame(term2_list, columns=["Terminator_Location"])
export_csv = term2_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/target_genome/term2_data.tsv", index=False)
term3_dataframe = pandas.DataFrame(term3_list, columns=["Terminator_Location"])
export_csv = term3_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/target_genome/term3_data.tsv", index=False)
rnase1_dataframe = pandas.DataFrame(data=rnase1_list, columns=["Rnase_Location"])
export_csv = rnase1_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/target_genome/rnase1_data.tsv", index=False)
rnase2_dataframe = pandas.DataFrame(data=rnase2_list, columns=["Rnase_Location"])
export_csv = rnase2_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/target_genome/rnase2_data.tsv", index=False)
rnase3_dataframe = pandas.DataFrame(rnase3_list, columns=["Rnase_Location"])
export_csv = rnase3_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/genome_coordinates/target_genome/rnase3_data.tsv", index=False)'''
