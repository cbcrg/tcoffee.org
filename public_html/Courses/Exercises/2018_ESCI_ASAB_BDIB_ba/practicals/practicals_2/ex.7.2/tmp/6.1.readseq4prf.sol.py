#!/usr/bin/env python

import sys
import re
class Seq:
    def __init__(self):
        self.id=0
        self.seq=""
        self.features=""
        
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

seq=readseq(sys.argv[1])
for name in (seq.keys()):
    sys.stdout.write(">%s\n%s\n"%(name, seq[name]))
