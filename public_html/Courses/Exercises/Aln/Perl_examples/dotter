#!/usr/bin/env perl

#
#   A course on sequence Alignments
# 	Cedric Notredame 2001
#   All rights reserved
#   Files can be redistributed without permission
#   Comercial usage is forbiden
#   ascii dotter


# Read The sequences from a fasta format file:
while (<>){$seq.=$_;}

#extract the names and the sequences
@name_list=($seq=~/>(.*)[^>]*/g);
@seq_list =($seq=~/>.*([^>]*)/g);

# get rid of the newlines, spaces and numbers
foreach $seq (@seq_list)
	{
	# get rid of the newlines, spaces and numbers
	$seq=~s/[\s\d]//g;	
	}
# split the sequences
for ($i=0; $i<=$#name_list; $i++)
	{
	$res[$i]=[$seq_list[$i]=~/([a-zA-Z]{1})/g];
	}

# make the dotter
print " ";
for ($j=0; $j<=$#{$res[1]}; $j++){print "$res[1][$j]";}
print "\n";

for ($i=0; $i<=$#{$res[0]}; $i++)
 	{
 	print "$res[0][$i]";
 	for ($j=0; $j<=$#{$res[1]}; $j++)
 		{
 		if ($res[1][$j] eq $res[0][$i]){print "*";$tot++;}
 		else {print " ";}
 		}
 	print "\n";
 	}
	
print "TOT=$tot\n";


                         
