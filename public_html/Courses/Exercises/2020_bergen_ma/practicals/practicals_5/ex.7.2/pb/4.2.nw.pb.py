#!/usr/bin/env python
import sys
import re
from Bio import SeqIO


def readmatrix(filename):
 #
 # This is where you want some matrix reading to happen
 #



#read the sequences
#read the sequences
record = list(SeqIO.parse(sys.argv[1], "fasta"))

#This is where you want to read your matrix

#remove gaps from sequences to align
seqI=re.sub(r'[\n-]',"",''.join(record[0].seq))
seqJ=re.sub(r'[\n-]',"",''.join(record[1].seq))

lenI=len(seqI)
lenJ=len(seqJ)

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
        if seqI[i-1]!='-' and seqJ[j-1]!='-':
            #This is where you want to get your score
        else:
            s=0

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


sys.stdout.write("Optimal Score: %d\n"%(int(smat[lenI][lenJ])))
i=lenI
j=lenJ
lenAln=0
alnI=[]
alnJ=[]

while ((i==0 and j==0)!=1):
    if (tb[i][j]==0):
        i-=1
        j-=1
        alnI.append(seqI[i])
        alnJ.append(seqJ[j])
    elif (tb[i][j]==-1):
        j-=1
        alnI.append('-')
        alnJ.append(seqJ[j])
    elif (tb[i][j]==1):
        i-=1
        alnI.append(seqI[i])
        alnJ.append('-')
    lenAln+=1

alnI=alnI[::-1]
alnJ=alnJ[::-1]

alnI=''.join(alnI)
alnJ=''.join(alnJ)

sys.stdout.write (">%s\n%s\n"%(record[0].id, alnI))
sys.stdout.write (">%s\n%s\n"%(record[1].id, alnJ))

