#!/usr/bin/env perl

#
#   A course on sequence Alignments
# 	Cedric Notredame 2001
#   All rights reserved
#   Files can be redistributed without permission
#   Comercial usage is forbiden
#   ktup dotter: compute the number of ktups for k=1


# Read The sequences from a fasta format file:
while (<>){$seq.=$_;}

#extract the names and the sequences
@name_list=($seq=~/>(.*)[^>]*/g);
@seq_list =($seq=~/>.*([^>]*)/g);

# get rid of the newlines, spaces and numbers
foreach $seq (@seq_list)
	{
	# get rid of the newlines, spaces numbers and gaps
	$seq=~s/[\s\d-]//g;	
	}
# split the sequences
for ($i=0; $i<=$#name_list; $i++)
	{
	$res[$i]=[$seq_list[$i]=~/([a-zA-Z]{1})/g];
	}

# Make the two ktup list
$alpha="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
@ktup_list=($alpha=~/([a-zA-Z]{1})/g);

for ($i=0; $i<=$#{$res[0]}; $i++)
  {
    $ktup_seq0{$res[0][$i]}{$nktup_seq0{$res[0][$i]}++}=$i;
  }
for ($j=0; $j<=$#{$res[1]}; $j++)
  {
    $ktup_seq0{$res[1][$j]}{$nktup_seq1{$res[1][$j]}++}=$j;
  }

# Enumerate matches
for ($a=0; $a<=$#ktup_list; $a++)
  {
    
    $ktup=$ktup_list[$a];
    
    $n=$nktup_seq0{$ktup};
    $m=$nktup_seq1{$ktup};
    $tot+=($n)*($m);
 
#   for ( $i=0; $i<=$n; $i++)
#	  {
#	    for ($j=0; $j<=$m; $j++)
#	      {
#		print "$ktup_seq0{$ktup}{$i}-$ktup_seq1{$ktup}{$j}}\n";
#	      }
#	  }
  }
print "TOT=$tot\n";
