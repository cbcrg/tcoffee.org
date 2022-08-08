#!/usr/bin/env perl

#   A course on sequence Alignments
# 	Cedric Notredame 2001
#   All rights reserved
#   Files can be redistributed without permission
#   Comercial usage is forbiden
#   dynamic programming profile-profile (special case, gop=0).
#   each sequence comes in a differrent file


#Parameters:
#matrix=idmat: 10 pour un Match, -10 pour un MM
$match=10;
$mismatch=-10;
$gop=-10;
$gep=-1;

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

$n0=$#name_list0+1;
$n1=$#name_list1+1;
$len0=$#{$res0[0]}+1;
$len1=$#{$res1[0]}+1;

#evaluate substitutions

for ($i=0; $i<=$len0; $i++){$smat[$i][0]=$i*$gep;$tb[$i][0 ]= 1;}
for ($j=0; $j<=$len1; $j++){$smat[0][$j]=$j*$gep;$tb[0 ][$j]=-1;}
	
for ($i=1; $i<=$len0; $i++)
	{
	for ($j=1; $j<=$len1; $j++)
		{
		#calcul du score
		for ($nm=0,$s=0,$n=0; $n<$n0; $n++)
		  {
		    for ($m=0; $m<$n1; $m++)
		      {
			if ( $res0[$n][$i-1] ne "-" && $res1[$m][$j-1] ne "-")
			  {
			    $nm++;
			    if ( $res0[$n][$i-1] eq $res1[$m][$j-1]){$s+=$match;}
			    else {$s+=$mismatch;}
			  }
		      }
		  }
		$s=$s/$nm;

		      
		$sub=$smat[$i-1][$j-1]+$s;
		$del=$smat[$i  ][$j-1]+$gep;
		$ins=$smat[$i-1][$j  ]+$gep;
		
		if   ($sub>$del && $sub>$ins){$smat[$i][$j]=$sub;$tb[$i][$j]=0;}
		elsif($del>$ins){$smat[$i][$j]=$del;$tb[$i][$j]=-1;}
		else {$smat[$i][$j]=$ins;$tb[$i][$j]=1;}
		}
	}

$i=$len0;
$j=$len1;
$aln_len=0;


while (!($i==0 && $j==0))
	{
	if ($tb[$i][$j]==0)
		{
		  $i--;
		  $j--;
		  for ( $n=0; $n<$n0; $n++){$aln0[$n][$aln_len]=$res0[$n][$i];}
		  for ( $m=0; $m<$n1; $m++){$aln1[$m][$aln_len]=$res1[$m][$j];}
		}
	elsif ($tb[$i][$j]==-1)
		{
		  $j--;
		  for ( $n=0; $n<$n0; $n++){$aln0[$n][$aln_len]='-';}
		  for ( $m=0; $m<$n1; $m++){$aln1[$m][$aln_len]=$res1[$m][$j];}
		}
	elsif ($tb[$i][$j]==1)
		{
		  $i--;
		  for ( $n=0; $n<$n0; $n++){$aln0[$n][$aln_len]=$res0[$n][$i];}
		  for ( $m=0; $m<$n1; $m++){$aln1[$m][$aln_len]='-';}
		}
	$aln_len++;	
	}
#Output en Fasta:
for ($n=0; $n<$n0; $n++)
  {
    print ">$name_list0[$n]\n";
    for ($i=$aln_len-1; $i>=0; $i--){print $aln0[$n][$i];}
    print "\n";
  }
for ( $m=0; $m<$n1; $m++)
  {
    print ">$name_list1[$m]\n";
    for ($j=$aln_len-1; $j>=0; $j--){print $aln1[$m][$j];}
    print "\n";
  }


                         
