import csv
import pandas


transcript_list = []
proteinY_list = []
proteinZ_list = []

for i in range(1, 151):
    transcript_list.append([i, "proteinX", (0.25 * i)])
    transcript_list.append([i, "proteinY", (0.3 * i)])
    if i == 150:
        proteinY_list = transcript_list[-1]
    transcript_list.append([i, "proteinZ", (0.1 * i)])
    if i == 150:
        proteinZ_list = transcript_list[-1]
for i in range(151, 241):
    transcript_list.append([i, "proteinX", (0.25 * i)])
    transcript_list.append([i, "proteinY", proteinY_list[-1]])
    transcript_list.append([i, "proteinZ", proteinZ_list[-1]])

transcript_dataframe = pandas.DataFrame()
export_csv = transcript_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/random_data.tsv", index=False)

file = open("random_data.tsv", "w")
writer = csv.writer(file, delimiter='\t')
writer.writerow(["time", "species", "transcript"])
for item in transcript_list:
    writer.writerow(item)
file.close()

print("Success!")
