import yaml

with open('gene.yml') as f:
    dataMap = yaml.safe_load(f)

dataMap['promoter1']['start'] = 0
print(dataMap['promoter1']['start'])
