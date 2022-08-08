#!/usr/sbin/perl
#
#   A course on sequence Alignments
# 	Cedric Notredame 2001
#   All rights reserved
#   Files can be redistributed without permission
#   Comercial usage is forbiden
#   score_aln


#Parameters: <aln>
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
$len_aln=$#{$res[1]};
$pstate=$score=0;

for ($i=0; $i<=$len_aln; $i++)
	{
	#Extension and substitution cost
	
	if    ($res[1][$i]eq '-' &&  $res[0][$i] eq '-'){;}
	elsif ($res[1][$i]eq '-'){$state=1;$score+=$gep}
	elsif ($res[0][$i]eq '-'){$state=2;$score+=$gep}
	else 
		{
		$state=3;
		if ($res[1][$i]eq $res[0][$i]){$score+=$match;}
		else {$score+=$mismatch;}
		}
	#Gop cost
	if      ($state==3){;}
	elsif ($state!=$pstate){$score+=$gop;}

	$pstate=$state;
	}	

print "SCORE=$score\n";	


                         
