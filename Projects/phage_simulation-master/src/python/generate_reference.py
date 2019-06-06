
'''
This script replaces the wild type gene 10A in the T7 reference genome with that of a recoded gene 10A provided in a fasta file.
'''

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

# Download T7 wild-type and E. coli genbank records
Entrez.email = "benjamin.r.jack@gmail.com"
handle = Entrez.efetch(db="nuccore",
                       id=["V01146.1", "U00096.3"],
                        rettype="gb",
                        retmode="text")

records = list(SeqIO.parse(handle, "genbank"))
phage_record = records[0]
ecoli_record = records[1]

# Pull out location of gene 10A
for feature in phage_record.features:
    if feature.type == 'CDS' and feature.qualifiers['protein_id'][0] == 'CAA24427.1':
        gene_location = feature.location

# Parse recoded gene
new_gene = SeqIO.read("recoded_10A.fasta", "fasta")
# Replace wildtype gene with recoded
new_sequence = phage_record.seq[0:gene_location.start] + new_gene.seq + phage_record.seq[gene_location.end:]
# Double check that the new sequence and old sequence are the same length
assert(len(phage_record.seq) == len(new_sequence))
# Generate new record for recoded T7
new_record = SeqRecord(new_sequence, 
                       id=phage_record.id, description=phage_record.description + new_gene.description)
# Write out wildtype and recoded records, with E. coli genome
recoded_records = [new_record, ecoli_record]
wildtype_records = [phage_record, ecoli_record]

SeqIO.write(recoded_records, "T7_recoded_ecoli_ref.fasta", "fasta")
SeqIO.write(wildtype_records, "T7_wildtype_ecoli_ref.fasta", "fasta")
