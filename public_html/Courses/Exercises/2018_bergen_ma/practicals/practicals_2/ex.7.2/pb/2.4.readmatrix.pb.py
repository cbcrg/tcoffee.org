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
  
  #get the alphabet and read it in the right order
  for ln in range (0,nlines-1):
      line=content[ln]
      splitted = line.split()
      a=splitted[0]
      if a not in matrix:
          matrix[a] = {}
          letters.append(a)
          
  #missing1 - Inspire yourself from the above block and get the values into matrix - about 10 lines         

         
  return (letters,matrix)


(alp,matrix)=readmatrix(sys.argv[1])

for x in range (0,len(alp)):
    aa1=alp[x]
    print "%c"%(aa1),
    for y in range (0,x+1):
        aa2=alp[y]
        print "%4d"%(int(matrix[aa1][aa2])),
    print "\n",

print "%-4c"%(' '),
for x in range (0,len(alp)):
    print "%-4c"%(alp[x]),

