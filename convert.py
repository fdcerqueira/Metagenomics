from Bio import SeqIO
import sys
import os


##parse from bash script
gbk_filename = sys.argv[1]
out_path = sys.argv[2]

#file with information extracted from .gbk file
path_file = os.sep.join([out_path,"snp.txt"])

#file stram to extract gene names, scaffold names,
#start position of gene in the scaffold, end start position of gene in the scaffold, strand, and the gene product

output = open(path_file, "w")

for seq_record in SeqIO.parse(gbk_filename, "genbank") :

    for seq_feature in seq_record.features :
        if seq_feature.type=="CDS":
            if 'gene' in seq_feature.qualifiers:

                line = "{};{};{};{};{};{}\n".format(seq_feature.qualifiers['gene'],
                                                  seq_record.name,
                                                  int(seq_feature.location.start),
                                                  int(seq_feature.location.end),
                                                  seq_feature.strand,seq_feature.qualifiers['product'])


                for ch in ["[","]","'"]:
                    line = line.replace(ch,"")

                output.write(line)

output.close()
