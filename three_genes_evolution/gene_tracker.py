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
writer.writerow(["geneX", 26, 121, 0, 0, 0])
writer.writerow(["geneY", 159, 280, 0, 0, 0])
writer.writerow(["geneZ", 319, 449, 0, 0, 0])
writer.writerow(["promoter1", 1, 10, 0.0, 0.0, 0])
writer.writerow(["promoter2", 0, 0, 0.0, 0.0, 0])
writer.writerow(["promoter3", 0, 0, 0.0, 0.0, 0])
writer.writerow(["rnase1", 0, 0, 0, 0, 0])
writer.writerow(["rnase2", 0, 0, 0, 0, 0])
writer.writerow(["rnase3", 0, 0, 0, 0, 0])
writer.writerow(["terminator1", 0, 0, 0.0, 0.0, 0])
writer.writerow(["terminator2", 0, 0, 0.0, 0.0, 0])
writer.writerow(["terminator3", 0, 0, 0.0, 0.0, 0])
writer.writerow(["f_old", 0, 0, 0, 0, 0.0])
writer.writerow(["f_new", 0, 0, 0, 0, 0.0])
writer.writerow(["region1", 11, 25, 0, 0, 0])
writer.writerow(["region2a", 122, 132, 0, 0, 0])
writer.writerow(["region2b", 133, 143, 0, 0, 0])
writer.writerow(["region2c", 144, 159, 0, 0, 0])
writer.writerow(["region3a", 281, 291, 0, 0, 0])
writer.writerow(["region3b", 292, 302, 0, 0, 0])
writer.writerow(["region3c", 303, 318, 0, 0, 0])

file.close()

print("Success!")
