#!/usr/bin/env python

import sys
import math
from Bio import SeqIO
from collections import defaultdict

records = list (SeqIO.parse(sys.argv[1], "fasta"))

alp="ACDEFGHIKLMNPQRSTVWY" # this is the AA alphabet


#initialize the hash with the full alphabet: you never know what's outhere
#set the default value to 1 to avoid overflow
D={}#will contain the total count/frequency of each transition and eventually the log odds
T={}#will contain the total count or frequency of each AA

#Start initilizing everything to 0 or to a pseudocount
totR=float (0);
totS=float (0);
for x in range (len(alp)):
    D[alp[x]]= {}
    T[alp[x]]=float(0)
    totR+=float(0)
    for y in range (len(alp)):
        D[alp[x]][alp[y]]=float(0)
        
#count all the AA occurences - mind the gaps!        
nseq=len(records)
for x in range (0,nseq):
    for z in range (len(records[x].seq)):
        aa1=records[x].seq[z].upper()
        if aa1 != '-':
            T[aa1]+=1
            totR+=1

#Turn the amino counts into relative frequencies 
for x in range (len(alp)):
    aa1=alp[x]
    T[aa1]/=totR

#count the number of transitions
for x in range (0,nseq-1):
    for y in range (x+1,nseq):
        for z in range (len(records[x].seq)):
            aa1=records[x].seq[z].upper()
            aa2=records[y].seq[z].upper()
            if aa1 != '-' and aa2 !='-':
                D[aa1][aa2]+=float (1)
                D[aa2][aa1]+=float (1)
                totS+=float (2)

# print the substitution matrix                
for x in range (0,len(alp)):
    aa1=alp[x]
    sys.stdout.write("%c"%(aa1))
    for y in range (0,x+1):
        aa2=alp[y]
        if T[aa1]>0 and T[aa2]>0 and D[aa1][aa2]>0:
            
            sub=(math.log10((D[aa1][aa2]/totS)/(T[aa1]*T[aa2])))*10
            sys.stdout.write("%4d"%(int(sub)))
        else:
            sys.stdout.write("%4d"%(0))
    sys.stdout.write("\n")

sys.stdout.write("%-4c"%(' '))
for x in range (0,len(alp)):
    sys.stdout.write("%-4c"%(alp[x]))
sys.stdout.write("\n")