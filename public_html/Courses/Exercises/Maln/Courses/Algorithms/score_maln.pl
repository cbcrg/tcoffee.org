#!/usr/bin/env perl
#
#   A course on sequence Alignments
# 	Cedric Notredame 2001
#   All rights reserved
#   Files can be redistributed without permission
#   Comercial usage is forbiden
#   score_aln


#Parameters:
#matrix=idmat: 10 pour un Match, -10 pour un MM
$match=10;
$mismatch=-10;
$gop=-10;
$gep=-1;

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
	$res[$i]=[$seq_list[$i]=~/([a-zA-Z-]{1})/g];
	}

#evaluate substitutions
$len_aln=$#{$res[1]}+1;
$nseq   =$#name_list+1;
$pstate=$score=0;

# Star Tree Evaluation
$alpha="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ-";
@symbol_list=($alpha=~/([a-zA-Z-]{1})/g);

for ($i=0; $i< $len_aln; $i++)
  {
    $best_col_score=$n;
    foreach $a (@symbol_list)
      {
	$col_score {$a}=0;	
	for ($n=0; $n<$nseq; $n++)
	  {	    
	    if ( $res[$n][$i]    eq "-" && $a eq "-"){;}
	    elsif ( $res[$n][$i] eq "-" || $a eq "-"){$col_score{$a}++;}
	    elsif ( $a eq $res[$n][$i]){;}
	    else {$col_score{$a}++;}
	  }
	if ( $col_score{$a}<$best_col_score){$best_col_score=$col_score{$a};}
      }
    $star_score+=$best_col_score;
  }

print "STAR SCORE=$star_score\n";

# DP Evaluation:

for ($i=0; $i< $len_aln; $i++)
  {
    for ($n=0; $n<$nseq-1; $n++)
      {
	for ($m=$n+1; $m<$nseq; $m++)
	  {
	    
	    if ( $res[$n][$i]    eq "-" && $res[$m][$i] eq "-"){;}
	    elsif ( $res[$n][$i] eq "-" || $res[$m][$i] eq "-"){$sp_score++;}
	    elsif ( $res[$m][$i] eq $res[$n][$i]){;}
	    else  { $sp_score++;}
	  }
      }
  }
print "SP   SCORE=$sp_score\n";
