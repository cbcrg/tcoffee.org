#!/usr/bin/env python
import sys
import re
from Bio import SeqIO

#read the sequences
record = list(SeqIO.parse(sys.argv[1], "fasta"))

#remove gaps from sequences to align
seqI=re.sub(r'[\n-]',"",''.join(record[0].seq))
seqJ=re.sub(r'[\n-]',"",''.join(record[1].seq))





lenI=len(seqI)
lenJ=len(seqJ)

match=1
mismatch=-10
gep=-4


smat = [[0 for x in range(lenJ+1)] for y in range(lenI+1)]
tb   = [[0 for x in range(lenJ+1)] for y in range(lenI+1)]

for i in range (0, lenI+1):
    smat[i][0]=i*gep
    tb[i][0]=1

for j in range (0, lenJ+1):
    smat[0][j]=j*gep
    tb[0][j]=-1


for i in range (1, lenI+1):
    for j in range (1, lenJ+1):
        if seqI[i-1]==seqJ[j-1]:
            s=match
        else:
            s=mismatch
        Sub=smat[i-1][j-1]+s
        Del=smat[i][j-1]+gep
        Ins=smat[i-1][j]+gep

        if Sub>Del and Sub >Ins:
            smat[i][j]=Sub
            tb  [i][j]=0  
        elif Del>Ins:
            smat[i][j]=Del
            tb[i][j]=-1
        else:
            smat[i][j]=Ins
            tb[i][j]=1

##missing1 - firgure out which cell contains the optimal score


