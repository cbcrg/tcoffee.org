#!/usr/bin/env python
import sys
import re

class Seq:
    def __init__(self):
        self.id=0
        self.seq=""
        self.features=""
        
#read lines one at a time
#do not store the whole file in memory if not needed
record=[]
nrec=-1
inseq=0
with open(sys.argv[1]) as f:
    for line in f:
        if re.match #missing1
            nrec+=1
            record.append(Seq())
            mobj=re.match #missing2

            if (mobj):
                record[nrec].id=mobj.group(1)
                record[nrec].features=mobj.group(2)
            inseq=0
        else :
            if inseq==0 :
                inseq=1
                #missing3
            else:
                cstring=record[nrec].seq+line
                #missing4
                
           
for x in range (0,nrec+1):
    print ">%s\n%s\n"%(record[x].id, record[x].seq),
