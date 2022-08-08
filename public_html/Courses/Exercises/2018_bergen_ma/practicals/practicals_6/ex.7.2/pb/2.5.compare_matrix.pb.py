#!/usr/bin/env python
import sys
import operator
from collections import defaultdict

#read the matrix
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


matrix1=readmatrix(sys.argv[1])
matrix2=readmatrix(sys.argv[2])

#get the alphabet
alp = []
for a in matrix1:
    alp.append(a)

#extrac the entries in   unsorted_mat dic and sort them in the sorted_mat ktuple list  
unsorted_mat={}
sorted_mat=()


for x in range (0,len(alp)):
    aa1=alp[x]
    for y in range (0,x+1):
        aa2=alp[y]
        delta=int(matrix1[aa1][aa2]) - int (matrix2[aa1][aa2])
        key=aa1+aa2
        unsorted_mat[key]=delta
        
sorted_mat=#missing1 - Figure out a way to sort the values
for x,y in sorted_mat:
  print (x, y)
