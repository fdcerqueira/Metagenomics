#!/usr/bin/python3

import itertools as it
import sys
import os

#parse input file, and outputpath from instrain.sh
input_file, output_path = sys.argv[1], sys.argv[2]

###join output path to final output file
path_file = os.sep.join([output_path,"compare_permutations.txt"])

#new list to remove \n"
new_list = []

#file stream input file
lines = open(input_file, "r")
content_list = lines.readlines()

#loop to remove \n
for i in content_list:
    new_list.append(i.strip())

print("")
print("Final list of samples")

lines.close()

#permutations between files without repetitions, to perform later inStrain compare
file=open(path_file, "w")
print("")
print("Sample's combinations:")

for a,b in it.combinations(new_list, 2):
    print("")
    print(a, b)
    file.write(a+" "+b+"\n")

file.close()
