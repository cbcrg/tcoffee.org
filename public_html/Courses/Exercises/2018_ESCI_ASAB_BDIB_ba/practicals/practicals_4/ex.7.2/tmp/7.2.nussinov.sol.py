#!/usr/bin/env python
import sys
import re

"""
The Nussinov Algorithm solves the problem of RNA non-crossing secondary
structure prediction by base pair maximization with input S.
"""

class Seq:
    def __init__(self):
        self.id=0
        self.seq=""
        self.features=""
                
def readseq (fname):
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
    for x in range (0,nrec+1):
	    record[x].seq=re.sub (r'[ \n\t\r]',"",record[x].seq)
    
    return record 


	    
def printStruc(seq,Struc):
	
	for k in range (len(seq)):
		sys.stdout.write (">%s\n%s\n"%(seq[k].id,seq[k].seq))
	sys.stdout.write (">nussinov\n%s\n"%(''.join(Struc)))
	
	

def multi_sigma (seq, i, j):
	n=len(seq)
	score=0
	
	for k in range (n):
		n1=seq[k].seq[i].upper()
		n2=seq[k].seq[j].upper()
		if (n1,n2) in (('A','U'),('U','A'),('C','G'),('G','C')):
			score+=1
		
	return score


def Nussinov(seq,N,L):
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
				N[i+1][j-1] + multi_sigma(seq,i,j),
				max([N[i][k] + N[k+1][j] for k in range(i,j)])
			]
                        N[i][j] = max(tmp)
        return N



def Traceback(seq,N,Struc,i,j):
	if j <= i: return (N, Struc)
	elif N[i][j] == N[i+1][j]:
		(N, Struc)=Traceback(seq,N,Struc,i+1, j)
	elif N[i][j]== N[i][j-1]:
		(N, Struc)=Traceback(seq,N,Struc,i,j-1)
	elif N[i][j]==(N[i+1][j-1]+multi_sigma(seq,i,j)):
		Struc[i] = '('
		Struc[j] = ')'
		(N, Struc)=Traceback(seq,N,Struc,i+1,j-1)
	else:
		for k in range(i+1,j-1):
			if N[i][j]==(N[i][k]+N[k+1][j]):
				(N, Struc)=Traceback(seq,N,Struc,k+1, j-1)
				(N, Struc)=Traceback(seq,N,Struc,i, k)
                                return (N, Struc)
	return (N, Struc)		



if __name__ == '__main__':
	seq=readseq(sys.argv[1])
	L=len(seq[0].seq)
	N=[ ['-']*L for i in range(L)]
	Struc = ['.']*L
	
	N=Nussinov(seq,N,L)
		
	(N, Struc)=Traceback(seq,N,Struc,0, L-1)
        
	printStruc(seq,Struc)
