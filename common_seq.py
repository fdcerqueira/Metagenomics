from Bio import SeqIO
import sys
import os

#parse input
common_fasta=sys.argv[1]
gbk_filename = sys.argv[2]
#parse output
output=sys.argv[3]


records = SeqIO.parse(common_fasta, "fasta")
idi = [record.id for record in records]

with open(output, "a") as f:
    final_records = []
    for record in SeqIO.parse(gbk_filename, "genbank"):
        for feature in record.features:
            if feature.type == "CDS" and feature.qualifiers["locus_tag"][0] in idi:
                pass
            else:
                record.features.remove(feature)
        SeqIO.write(record, f, "genbank")
