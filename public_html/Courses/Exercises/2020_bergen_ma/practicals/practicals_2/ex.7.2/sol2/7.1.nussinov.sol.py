#!/usr/bin/env python
import sys
"""
The Nussinov Algorithm solves the problem of RNA non-crossing secondary
structure prediction by base pair maximization with input S.
"""

Seq = 'cccccccggccggcauggucccagccuccucgcuggcgccggcugggcaacaccauugcacuccgguggcgaaugggggg'

L = len(Seq)
N = [ ['-']*L for i in range(L)]
Struc = ['.']*L



def printStruc():
	str=""
	sys.stdout.write ("'[Secondary Structure]'\n")
	sys.stdout.write("%c%s\n"%('-',Seq))
	sys.stdout.write("%c%s"%('-',str.join(Struc)))

def sigma(i,j):
	n1=Seq[i].upper()
	n2=Seq[j].upper()
	if (n1, n2) in (('A','U'),('U','A'),('C','G'),('G','C')):
		return 1
	return 0

def Nussinov():
	# initialization
	for i in range(L):
		N[i][i] = 0
		if i != 0: N[i][i-1] = 0
	# recursion
	for a in range(1,L):
		for i in range(L-a):
			j = i + a
			tmp = [
				N[i+1][j],
				N[i][j-1],
				N[i+1][j-1] + sigma(i,j),
				max([N[i][k] + N[k+1][j] for k in range(i,j)])
			]
			N[i][j] = max(tmp)
	print ("SCORE=",N[0][L-1])

def Traceback(i,j):
	if j <= i: return
	elif N[i][j] == N[i+1][j]:
		Traceback(i+1, j)
	elif N[i][j]== N[i][j-1]:
		Traceback(i,j-1)
	elif N[i][j]==(N[i+1][j-1]+sigma(i,j)):
		Struc[i] = '('
		Struc[j] = ')'
		Traceback(i+1,j-1)
	else:
		for k in range(i+1,j-1):
			if N[i][j]==(N[i][k]+N[k+1][j]):
				Traceback(k+1, j-1)
				Traceback(i, k)
				return
	return		
if __name__ == '__main__':
	Nussinov()
	Traceback(0, L-1)
	printStruc()
