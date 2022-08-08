#!/usr/bin/env python

import sys
import re
class Seq:
    def __init__(self):
        self.id=0
        self.seq=""
        self.features=""



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

def align (seq,Igroup, Jgroup, matrix, gep):
    
    lenI=#get then length of the first group
    lenJ=#get the length of the second group
    
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
        s=0
        nsub=0
        #missing 1 - get s, the average score af all against all matches between the two columns
                         
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


    sys.stdout.write ("Optimal Score: %d\n"%(int(smat[lenI][lenJ])))

         

def readseq (fname):
    seqlist={}
    record=[]
    nrec=-1
    inseq=0
    with open(fname) as f:
      for line in f:
        if re.match ( r'^>', line):
            nrec+=1
            record.append(Seq())
            mobj=re.match ( r'^>(\S*)\s*(.*)', line)

            if (mobj):
                record[nrec].id=mobj.group(1)
                record[nrec].features=mobj.group(2)
            inseq=0
        else :
            if inseq==0 :
                inseq=1
                record[nrec].seq=line
            else:
                cstring=record[nrec].seq+line
                record[nrec].seq=cstring
                
    seqlist={}       
    for x in range (0,nrec+1):
        record[x].seq=re.sub (r'[ \n\t\r]',"",record[x].seq)
        seqlist[record[x].id]=record[x].seq

    return seqlist 
     



#The main code starts here

gep=-4
left=readseq(sys.argv[1])
right=readseq(sys.argv[2])
matrix=readmatrix (sys.argv[3])
seq={}
seq=left.copy()
seq.update(right)
align (seq,list(left.keys()), list(right.keys()), matrix, gep)
