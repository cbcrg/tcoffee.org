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

#missing1 - Start initializing everything to 0 or to a pseudocount - can be done in 7-8 lines

        
#missing2 - count all the AA occurences - mind the gaps! - 6 lines        

#Turn the amino counts into relative frequencies 
for x in range (len(alp)):
    aa1=alp[x]
    T[aa1]/=totR

#missing3 - count the number of transitions - these are the importaat counts. Check in the next section to see how these values arre being handled when the log-odd is computed - 9 line


#Compute and Print the substitution matrix                
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