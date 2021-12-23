import pandas as pd
import sys
import os
import warnings
warnings.filterwarnings('ignore')


    #input
fich = sys.argv[1]
table = sys.argv[2]
tax = sys.argv[3]
meta_bin = sys.argv[4]
stb = sys.argv[5]
    #name for final table
name_final_table = sys.argv[6]
    #working directory
out_f=sys.argv[7]


path_file = os.sep.join([out_f,os.path.basename(name_final_table)])

    ##stb
stb=pd.read_table(stb, sep="\t", header= None)

stb.columns= ["scaffold","bin"]
stb['bin'] = stb['bin'].str.replace(r'.rename.fa', '')

    ##taxanomy database
fich=pd.read_table(fich, sep="\t")

fich=fich.rename(columns = {'user_genome':'bin'})
fich=fich[['bin', 'classification']]

    # new data frame with split value columns
fich1=fich['classification'].str.split(';', 7, expand=True)
fich1.rename(columns={0:"Kingdom",1:"Phylum",2:"Class",3:"Order",4:"Family",5:"Genus",6:"Species"}, inplace=True)

    #fucntion to remove the character "*__"
def rem(x):
    return x.str[3:]

    ##apply function
fich2 = fich1.apply(rem, axis=1)

    ##drop classification column from dataframe
fich=fich.drop(["classification"], axis=1)

    #cbind fich1 to bins
tax_database = pd.concat([fich2.reset_index(drop=True),fich.reset_index(drop=True)], axis=1)

    #merge stb with database taxonomy
al=tax_database.merge(stb, on='bin', how="right")

    #bins taxonomy
tax=pd.read_table(tax, sep="\t", header=None)
tax=tax.rename(columns={0:"scaffold",1:"bin"})

    ##keep numbers of scaffolds
tax1=tax['scaffold'].str.extract('(\d+)').astype(int)

    #add new column and remeve old scaffold
tax=tax.drop(["scaffold"], axis=1)
taxa_new = pd.concat([tax.reset_index(drop=True),tax1.reset_index(drop=True)], axis=1)
taxa_new=taxa_new.rename(columns={0:"scaffold"})

    #meta_bin
meta_bin=pd.read_table(meta_bin, sep="\t")

    #function to remove the two type of extensions
def remo(x):
    return x.str.replace('.contigs',"")

#apply the function
meta_bin2 = meta_bin.apply(remo, axis=1)
meta_bin2=meta_bin2.rename(columns={"Unnamed: 0":"bin", "superkingdom":"Kingdom","phylum":"Phylum","class":"Class","order":"Order","family":"Family","genus":"Genus","species":"Species"})

    #merge bins/scaffolds with taxonomy
al_bin=taxa_new.merge(meta_bin2, on='bin', how="right")

    #change position of columns
al_bin = al_bin[["Kingdom","Phylum","Class","Order","Family","Genus","Species","bin","scaffold"]]

    #append database+bins
all_t=al_bin.append(al)
    #merge taxonomy with final table of snps
table=pd.read_csv(table)
final_table=table.merge(all_t, on="scaffold", how="left")
#save final table
final_table.to_csv(path_file, index=False, na_rep="NA")
