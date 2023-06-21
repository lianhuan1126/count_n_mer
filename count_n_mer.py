import os
import sys
import re
from Bio import SeqIO
import csv
from collections import Counter

#%%
root_dir = str(sys.argv[1])
sys.path.append(root_dir)

#%% read files
file_path = root_dir+str(sys.argv[2])
out_path = root_dir+str(sys.argv[3])
mer_count =root_dir+str(sys.argv[4])


mer_count = int(sys.argv[4])

# read the sequences from the FASTA file
sequences = []
with open(file_path, "r") as f:
    for record in SeqIO.parse(f, "fasta"):
        sequences.append(str(record.seq))

# count all possible 6-mers in the sequences
kmers = Counter()
for seq in sequences:
    for i in range(len(seq) - mer_count + 1):
        kmer = seq[i:i+mer_count]
        kmers[kmer] += 1

# write the kmer counts to a file

with open(out_path, "w") as f:
    for kmer, count in kmers.items():
        f.write(f"{kmer}\t{count}\n")
