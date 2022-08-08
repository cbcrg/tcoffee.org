#!/usr/bin/env python
import sys
import math
from Bio import SeqIO
from collections import defaultdict

records = list (SeqIO.parse(sys.argv[1], "fasta"))
Nseq=len(records)
Laln=len (records[0].seq)

gcount=[0]*Laln
for c in range (0, Laln):
    for s in range (0, Nseq):
        letter=records[s].seq[c]
        if (letter == '-'):
            gcount[c]+=1;

for c in range (0, Laln):
    gcount[c]/=float(Nseq)
    gcount[c]*=100

for s in range (0, Nseq):
    print ">%s\n"%(records[s].id),
    for c in range (0, Laln):
        if gcount[c]<=int(sys.argv[2]):
            sys.stdout.write("%c"%(records[s].seq[c]))
    sys.stdout.write("\n")                        
        
    
