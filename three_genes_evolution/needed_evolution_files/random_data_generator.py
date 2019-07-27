import csv
import pandas


transcript_list = []
proteinX_list = []
proteinY_list = []
proteinZ_list = []

for i in range(1, 151):
    transcript_list.append([i, "proteinX", (0.2 * i)])
    if i == 150:
        proteinX_list = transcript_list[-1]
    transcript_list.append([i, "proteinY", (0.3 * i)])
    if i == 150:
        proteinY_list = transcript_list[-1]
    transcript_list.append([i, "proteinZ", (0.4 * i)])
    if i == 150:
        proteinZ_list = transcript_list[-1]
for i in range(151, 241):
    transcript_list.append([i, "proteinX", (0.2 * i)])
    transcript_list.append([i, "proteinY", (0.3 * i)])
    transcript_list.append([i, "proteinZ", (0.4 * i)])

transcript_dataframe = pandas.DataFrame()
export_csv = transcript_dataframe.to_csv("~/pinetree-toys/three_genes_evolution/needed_evolution_files/random_data3.tsv", index=False)

file = open("random_data3.tsv", "w")
writer = csv.writer(file, delimiter='\t')
writer.writerow(["time", "species", "transcript"])
for item in transcript_list:
    writer.writerow(item)
file.close()

print("Success!")
