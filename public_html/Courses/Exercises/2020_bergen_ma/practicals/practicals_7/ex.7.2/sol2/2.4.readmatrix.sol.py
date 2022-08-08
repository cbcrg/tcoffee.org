#!/usr/bin/env python
import sys
from collections import defaultdict

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
         
  return (letters,matrix)


(alp,matrix)=readmatrix(sys.argv[1])

for x in range (0,len(alp)):
    aa1=alp[x]
    sys.stdout.write ("%c"%(aa1))
    for y in range (0,x+1):
        aa2=alp[y]
        sys.stdout.write("%4d"%(int(matrix[aa1][aa2])))
    sys.stdout.write ("\n")

sys.stdout.write("%-4c"%(' '))
for x in range (0,len(alp)):
    sys.stdout.write("%-4c"%(alp[x]))
