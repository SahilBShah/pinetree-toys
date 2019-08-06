import yaml

with open('new_gene.yml') as f:
    genome_tracker = yaml.safe_load(f)

print(genome_tracker['rnase1']['start'])
