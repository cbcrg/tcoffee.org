#!/usr/bin/env python
import sys
import re

def readmatrix(filename):
  handle = open(filename, "r")
  content = handle.readlines()
  handle.close()

  matrix = {}
  letters = []
  nlines=len (content)
  
  
  for ln in range (0,nlines-1):
      line=content[ln]
      splitted = line.split()
      a=splitted[0]
      if a not in matrix:
          matrix[a] = {}
          letters.append(a)
          
  for ln in range (0,nlines-1):
      line=content[ln]
      splitted = line.split()
      l=len(splitted)
      aa1=splitted[0]
      for a in range (1,l):
          aa2=letters[a-1]
          matrix[aa1][aa2]=splitted[a]
          matrix[aa2][aa1]=splitted[a]
         
  return matrix




#read the sequences
record = list(SeqIO.parse(sys.argv[1], "fasta"))
matrix = readmatrix (sys.argv[2])

#remove gaps from sequences to align
seqI=re.sub(r'[\n-]',"",''.join(record[0].seq))
seqJ=re.sub(r'[\n-]',"",''.join(record[1].seq))

lenI=len(seqI)
lenJ=len(seqJ)

gep=-4


smat = [[0 for x in range(lenJ+1)] for y in range(lenI+1)]
tb   = [[0 for x in range(lenJ+1)] for y in range(lenI+1)]

bscore=0
for i in range (0, lenI+1):
    smat[i][0]=## missing 1- Can a local aln start with a gap?
    tb[i][0]=## missing 2 
    
      
for j in range (0, lenJ+1):
    smat[0][j]=## missing 3
    tb[0][j]= ## missing 4 

bscore=0
bi=0
bj=0

for i in range (1, lenI+1):
    for j in range (1, lenJ+1):
        if seqI[i-1]!='-' and seqJ[j-1]!='-':
            s=int(matrix[seqI[i-1]][seqJ[j-1]])
        else:
            s=0

        Sub=smat[i-1][j-1]+s
        Del=smat[i][j-1]+gep
        Ins=smat[i-1][j]+gep

        if Sub>Del and Sub >Ins ## missing 5
          smat[i][j]=Sub
          tb  [i][j]=0  
        elif Del>Ins ## missing 6
          smat[i][j]=Del
          tb[i][j]=-1
        #missing 7
          smat[i][j]=Ins
          tb[i][j]=1
        ##missing 8

        ##missing 9 - you need to remeber where to start the traceback from

          
print "Optimal Score: %d\n"%(int(smat[lenI][lenJ]))
i=##where do we start
j=##where do we start
lenAln=0
alnI=[]
alnJ=[]

while (tb[i][j] ## missing 10 - whwn will the traceback be finished?):
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

