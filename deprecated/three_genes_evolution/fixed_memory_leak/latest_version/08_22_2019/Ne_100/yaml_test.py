import yaml
import save

with open('new_gene.yml') as f:
    genome_tracker = yaml.safe_load(f)

genome_tracker['rnase1']['start'] = 5
print(genome_tracker['rnase1']['start'])
with open('new_gene.yml', 'w') as f:
    yaml.dump(genome_tracker, f)

f.close()
