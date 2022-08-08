#!/usr/bin/env perl

#
#   A course on sequence Alignments
# 	Cedric Notredame 2001
#   All rights reserved
#   Files can be redistributed without permission
#   Comercial usage is forbiden
#   dna2rna


# Read The sequences from a fasta format file:
# $seq=">name1\nATGCTA\nGATC\n>name2\nAGAC\nTCAT\n";

while (<>){$seq=$seq.$_;}

#extract the names and the sequences
@name_list=($seq=~/>(.*)[^>]*/g);
@seq_list =($seq=~/>.*([^>]*)/g);


foreach $seq (@seq_list)
	{
	# get rid of the newlines, spaces and numbers
	$seq=~s/[\s\d]//g;
	#Traduction
	$seq=~tr/AGCTagct/UCGAucga/;		
	}
#Print the result in Fasta Format

for ($i=0; $i<=$#name_list; $i++)
	{
	print ">$name_list[$i]\n$seq_list[$i]\n";
	}

                         
