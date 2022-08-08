#!/usr/bin/env python

import sys
import re
class Seq:
    def __init__(self):
        self.id=""
    def __init__(self):
        self.seq=""
    def __init__(self):
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
                
    # missing1 - add all the sequences in a dictionnary here - 3 lines
    return seqlist

#The main code starts here

seq=readseq(sys.argv[1])
for name in (seq.keys()):
    sys.stdout.write(">%s\n%s\n"%(name, seq[name]))
