#!/usr/bin/env python
import sys
from Bio import SeqIO

for record in SeqIO.parse(sys.argv[1], "fasta"):
    sys.stdout.write("%s\n"%(record.seq))

