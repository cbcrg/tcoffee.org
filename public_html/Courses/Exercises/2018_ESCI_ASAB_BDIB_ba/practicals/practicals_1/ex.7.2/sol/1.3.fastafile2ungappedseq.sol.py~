#!/usr/bin/env python
import sys
import re

class Seq:
    def __init__(self):
        self.id="x"
    def __init__(self):
        self.seq="y"
    def __init__(self):
        self.features="z"
        

#read lines one at a time
#do not store the whole file in memory if not needed
record=[]
nrec=-1
inseq=0
with open(sys.argv[1]) as f:
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
    string=re.sub(r'[\n-]',"",record[x].seq)
    record[x].seq=string
    print ">%s\n%s\n"%(record[x].id, record[x].seq),
