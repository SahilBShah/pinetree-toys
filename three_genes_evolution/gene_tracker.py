import csv
import pandas as pd

empty_list1 = []
empty_list2 = []
empty_list3 = []
empty_list4 = []
empty_list5 = []

tracker_dataframe = pd.DataFrame()
export_csv = tracker_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/gene_tracker.tsv", index=False)

file = open("gene_tracker.tsv", "w")
writer = csv.writer(file, delimiter='\t')
writer.writerow(["species", "start", "stop", "new_strength", "previous_strength", "value"])
writer.writerow(["geneX", 26, 121, "N/A", "N/A", "N/A"])
writer.writerow(["geneY", 159, 280, "N/A", "N/A", "N/A"])
writer.writerow(["geneZ", 319, 449, "N/A", "N/A", "N/A"])
writer.writerow(["promoter1", 1, 10, 0.0, 0.0, "N/A"])
writer.writerow(["promoter2", 0, 0, 0.0, 0.0, "N/A"])
writer.writerow(["promoter3", 0, 0, 0.0, 0.0, "N/A"])
writer.writerow(["rnase1", 0, 0, "N/A", "N/A", "N/A"])
writer.writerow(["rnase2", 0, 0, "N/A", "N/A", "N/A"])
writer.writerow(["rnase3", 0, 0, "N/A", "N/A", "N/A"])
writer.writerow(["terminator1", 0, 0, 0.0, 0.0, "N/A"])
writer.writerow(["terminator2", 0, 0, 0.0, 0.0, "N/A"])
writer.writerow(["terminator3", 0, 0, 0.0, 0.0, "N/A"])
writer.writerow(["f_old", "N/A", "N/A", "N/A", "N/A", 0.0])
writer.writerow(["f_new", "N/A", "N/A", "N/A", "N/A", 0.0])
writer.writerow(["region1", 11, 25, "N/A", "N/A", "N/A"])
writer.writerow(["region2a", 122, 132, "N/A", "N/A", "N/A"])
writer.writerow(["region2b", 133, 143, "N/A", "N/A", "N/A"])
writer.writerow(["region2c", 144, 159, "N/A", "N/A", "N/A"])
writer.writerow(["region3a", 281, 291, "N/A", "N/A", "N/A"])
writer.writerow(["region3b", 292, 302, "N/A", "N/A", "N/A"])
writer.writerow(["region3c", 303, 318, "N/A", "N/A", "N/A"])

file.close()

print("Success!")
