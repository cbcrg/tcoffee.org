#!/usr/bin/env python
import sys
import math
from Bio import SeqIO
from Bio import AlignIO
from Bio import Alphabet
from Bio.Alphabet import IUPAC
from Bio.Align import AlignInfo
from collections import defaultdict
from Bio import SubsMat

#This script will compute a log-odd matrix using bio-python



filename=sys.argv[1]
alpha = Alphabet.Gapped(IUPAC.protein)
c_align = AlignIO.read(filename, "clustal", alphabet=alpha)
summary_align = AlignInfo.SummaryInfo(c_align)

m0=summary_align.replacement_dictionary()
m1=SubsMat.SeqMat(m0)
m2=SubsMat.make_log_odds_matrix(m1,logbase=10,factor=10)
m2.print_mat()
