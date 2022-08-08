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
$gep=-10;
@BLOSUM62 = (
	[4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0,-2,-1,0,-4],
	[-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3,-1,0,-1,-4],
	[-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3,3,0,-1,-4],
	[-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3,4,1,-1,-4],
	[0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2,-4],
	[-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2,0,3,-1,-4],
	[-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2,1,4,-1,-4],
	[0,-2,0,-1,-3,-2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3,-1,-2,-1,-4],
	[-2,0,1,-1,-3,0,0,-2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3,0,0,-1,-4],
	[-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3,-1,3,-3,-3,-1,-4],
	[-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1,-4,-3,-1,-4],
	[-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2,0,1,-1,-4],
	[-1,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1,-3,-1,-1,-4],
	[-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1,-3,-3,-1,-4],
	[-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2,-2,-1,-2,-4],
	[1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2,0,0,0,-4],
	[0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0,-1,-1,0,-4],
	[-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3,-4,-3,-2,-4],
	[-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1,-3,-2,-1,-4],
	[0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4,-3,-2,-1,-4],
	[-2,-1,3,4,-3,0,1,-1,0,-3,-4,0,-3,-3,-2,0,-1,-4,-3,-3,4,1,-1,-4],
	[-1,0,0,1,-3,3,4,-2,0,-3,-3,1,-1,-3,-1,0,-1,-3,-2,-2,1,4,-1,-4],
	[0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2,0,0,-2,-1,-1,-1,-1,-1,-4],
	[-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,1]
	     );
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
		    #$s= getSubstScore(\@BLOSUM62,$res0[0][$i-1],$res1[0][$j-1]);  
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

print "$best_score\n";

#Traceback Procedure
$i=$best_i;
$j=$best_j;
$aln_len=0;
while (!($i==0 && $j==0) && $tb[$i][$j]!=2)
	{
	if ($tb[$i][$j]==0)
		{
		$aln0[$aln_len]=$res0[0][--$i];
		$aln1[$aln_len]=$res1[0][--$j];
		}
	elsif ($tb[$i][$j]==-1)
		{
		$aln0[$aln_len]='-';
		$aln1[$aln_len]=$res1[0][--$j];
		}
	elsif ($tb[$i][$j]==1)
		{
		$aln0[$aln_len]=$res0[0][--$i];
		$aln1[$aln_len]='-';
		}
	$aln_len++;
	
	}
$i++;
$j++;
#Output en Fasta:
print ">$name_list0[0]: [$i-$best_i]\n";
for ($i=$aln_len-1; $i>=0; $i--){print $aln0[$i];}
print "\n";
print ">$name_list1[0]: [$j-$best_j]\n";
for ($j=$aln_len-1; $j>=0; $j--){print $aln1[$j];}
print "\n";



                         


sub getSubstScore
{
	$aaOrder = "ARNDCQEGHILKMFPSTWYV";
	$matrix = shift;
	$aa1 = shift;
	$aa2 = shift;

	$pos1 = index($aaOrder, $aa1);
	$pos2 = index($aaOrder, $aa2);

	return ( $matrix->[$pos1][$pos2] );
}



