import pandas as pd
import sys
import os
import warnings
warnings.filterwarnings('ignore')
 #input
table = sys.argv[1]
meta_bin = sys.argv[2]
stb = sys.argv[3]
    #name for final table
name_final_table = sys.argv[4]
    #working directory
out_f=sys.argv[5]


path_file = os.sep.join([out_f,os.path.basename(name_final_table)])

stb=pd.read_table(stb, sep="\t", header= None)

stb.columns= ["scaffold","bin"]
stb['bin'] = stb['bin'].str.replace(r'.rename.fa', '')

#meta_bin
meta_bin=pd.read_table(meta_bin, sep="\t")

#function to remove the two type of extensions
def remo(x):
    return x.str.replace('.contigs',"")

#apply the function
meta_bin2 = meta_bin.apply(remo, axis=1)
meta_bin2=meta_bin2.rename(columns={"Unnamed: 0":"bin", "superkingdom":"Kingdom","phylum":"Phylum","class":"Class","order":"Order","family":"Family","genus":"Genus","species":"Species"})

al_bin=meta_bin2.merge(stb, on='bin', how="right")

#change position of columns
al_bin = al_bin[["Kingdom","Phylum","Class","Order","Family","Genus","Species","bin","scaffold"]]

#merge taxonomy with final table of snps

table=pd.read_csv(table)
final_table=table.merge(al_bin, on='scaffold', how="left")

#save final table
final_table.to_csv(path_file, index=False, na_rep="NA")
