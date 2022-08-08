#!/usr/bin/env python


from Bio import SeqIO

for record in SeqIO.parse(sys.argv[0], "fasta"):
    print(record.seq)

