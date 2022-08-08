#!/usr/bin/env perl

#   A course on sequence Alignments
# 	Cedric Notredame 2001
#   All rights reserved
#   Files can be redistributed without permission
#   Comercial usage is forbiden
#   Progressive Alignment: Sequences are aligned one by one according to the incoming order.
#   This programm calls the prf-prf alignment routine pnw.pl

#Parameters: <sequences>

# Read The sequences from a fasta format file:
while ( <>){$al.=$_;}

#extract the names and the sequences
@name_list=($al=~/>(.+)[^>]+/g);
@seq_list =($al=~/>.+([^>]+)/g);
$nseq=$#name_list+1;

#prepare sequence files to be sent to dpp.pl
open TMP1, ">tmp1";
print TMP1 ">$name_list[0]\n$seq_list[0]";
close (TMP1);

for ( $a=1; $a< $nseq; $a++)
  {
    open TMP2, ">tmp2";
    print TMP2 ">$name_list[$a]\n$seq_list[$a]";
    close (TMP2);

    $sub_aln=`pnw.pl tmp1 tmp2>tmp3`;
    `mv tmp3 tmp1`;
  }

open TMP, "tmp1";
while (<TMP>){print "$_";}

unlink ("tmp1");
unlink ("tmp2");
    

