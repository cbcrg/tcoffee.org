#!/usr/bin/env perl

#   A course on sequence Alignments
# 	Cedric Notredame 2001
#   All rights reserved
#   Files can be redistributed without permission
#   Comercial usage is forbiden
#   Smith et Waterman (special case, gop=0).
#   each sequence comes in a differrent file


#Parameters:
#matrix=idmat: 10 pour un Match, -10 pour un MM
$match=10;
$mismatch=-10;
$gop=-10;

# Read The two sequences from two fasta format file:



open F0, $ARGV[0];
while ( <F0>){$al0.=$_;}
close F0;

open F1, $ARGV[1];
while ( <F1>){$al1.=$_;}
close F1;


#extract the names and the sequences
@name_list0=($al0=~/>(.*)[^>]*/g);
@seq_list0 =($al0=~/>.*([^>]*)/g);

@name_list1=($al1=~/>(.*)[^>]*/g);
@seq_list1 =($al1=~/>.*([^>]*)/g);

# get rid of the newlines, spaces and numbers
foreach $seq (@seq_list0)
	{
	# get rid of the newlines, spaces and numbers
	$seq=~s/[\s\d]//g;	
	}
foreach $seq (@seq_list1)
	{
	# get rid of the newlines, spaces and numbers
	$seq=~s/[\s\d]//g;	
	}

# split the sequences
for ($i=0; $i<=$#name_list0; $i++)
	{
	$res0[$i]=[$seq_list0[$i]=~/([a-zA-Z-]{1})/g];
	}
for ($i=0; $i<=$#name_list1; $i++)
	{
	$res1[$i]=[$seq_list1[$i]=~/([a-zA-Z-]{1})/g];
	}

print ">$name_list0[0]\n$seq_list1[0]\n";
print ">$name_list0[0]\n$seq_list1[0]\n";
