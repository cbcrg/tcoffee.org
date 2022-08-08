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

#Smith and Waterman Recursion: Find the score of the optimal alignment
$len0=$#{$res0[0]}+1;
$len1=$#{$res1[0]}+1;
$zero=0;
for ($i=0; $i<=$len0; $i++){$smat[$i][0]=$zero;$tb[$i][0 ]=2;}
for ($j=0; $j<=$len1; $j++){$smat[0][$j]=$zero;$tb[0 ][$j]=2;}
	
for ($i=1; $i<=$len0; $i++)
	{
	for ($j=1; $j<=$len1; $j++)
		{
		#calcul du score
		if ($res0[0][$i-1] eq $res1[0][$j-1]){$s=$match;}
		else {$s=$mismatch;}
		
		$sub=$smat[$i-1][$j-1]+$s;
		$del=$smat[$i  ][$j-1]+$gep;
		$ins=$smat[$i-1][$j  ]+$gep;
		
		if   ($sub>$del && $sub>$ins && $sub> $zero){$smat[$i][$j]=$sub;$tb[$i][$j]=0;}
		elsif($del>$ins && $del>$zero              ){$smat[$i][$j]=$del;$tb[$i][$j]=-1;}
		elsif($ins>$zero                           ){$smat[$i][$j]=$ins;$tb[$i][$j]=1;}
		else {$smat[$i][$j]=$zero;$tb[$i][$j]=2;}
	
		if ($smat[$i][$j]> $best_score)
		  {
		     $best_score=$smat[$i][$j];
		     $best_i=$i;
		     $best_j=$j;
		   }
	        }
    }

print "BEST SCORE=$best_score\n";
