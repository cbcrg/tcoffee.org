#!/usr/bin/env python
import sys
import math
from Bio import SeqIO
from collections import defaultdict

records = list (SeqIO.parse(sys.argv[1], "fasta"))
Nseq=len(records)
Laln=len (records[0].seq)

gcount=[0]*Laln
#missing1 - count all the gaps on every position - can be done in 5 lines


for c in range (0, Laln):
    gcount[c]/=float(Nseq)
    gcount[c]*=100

#missing 2 - print the sequences one column at a time taking count into consideration - can be done in 6 lines
        
    
