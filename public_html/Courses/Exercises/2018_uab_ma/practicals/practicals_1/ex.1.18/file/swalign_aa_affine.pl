#!/usr/bin/perl


#
# This script reads two fasta-formatted sequence files and performs a basic version of the
# Smith-Waterman algorithm following approximately the recipe given in Durbin et al. 1998
# "Biological Sequence Analysis - Probabilistic models of proteins and nucleic acids",
# pp. 22-24.
#

# List of global variables:
# @x  seq 1
# $XL length of 1st sequence
# @y  seq 2
# $YL length of 2nd sequence
# @s  scoring matrix
# @M  match alignment matrix
# @I  insert alignment matrix

if (scalar(@ARGV) != 4) {die "Usage: ./script-name <sequence file 1> <sequence file 2> <gop> <gep>" ;}

# Read sequences and gap penalties from the command line.
$fasta1 = @ARGV[0];
$fasta2 = @ARGV[1];
$gop	= @ARGV[2]; # gop = Gap Opening Penalty
$gep	= @ARGV[3]; # gep = Gap Extension Penalty

# Substitution Matrix initialization (Blosum62 from the EMBOSS package)
#  Matrix made by matblas from blosum62.iij
#  * column uses minimum score
#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
#  Blocks Database = /data/blocks_5.0/blocks.dat
#  Cluster Percentage: >= 62
#  Entropy =   0.6979, Expected =  -0.5209
#A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
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




# input parsing
$/ = "\n>";                                # entry delimiter for Fasta-formatted sequence library
open(SEQ1, "$fasta1");                    # reading the first sequence
while(<SEQ1>)                              # major loop over first sequence in the input file
{
  /^>?([^\n]*)\n(.*)/s;                    # header goes into $1, sequence into $2
  $name = $1;
  $seq = $2;                               # otherwise $2 would be erased during next pattern match
  $seq =~ s/[^ARNDCQEGHILKMFPSTWYV]//gs;   # delete all characters in $seq except amino acids: ARNDCQEGHILKMFPSTWYV
  @x = split //, " " . $seq;               # make sequence array, note $x[0] = ' ', $x[1] = 1st base
  $XL = scalar(@x) - 1;                    # $XL is lengt of first sequence
  print "\nSequence 1: $name\n\n";         # print header of first sequence
  print "  ", @x, "\n";                    # print first sequence
}
close(SEQ1);

open(SEQ2, "$fasta2");                    # reading the first sequence
while(<SEQ2>)                              # same as for the first sequence
{
  /^>?([^\n]*)\n(.*)/s;
  $name = $1;
  $seq = $2;
  $seq =~ s/[^ARNDCQEGHILKMFPSTWYV]//gs;
  @y = split //, " " . $seq;
  $YL = scalar(@y) - 1;                    # $YL is lengt of second sequence
  print "\nSequence 2: $name\n\n";         # print header of second sequence
  print "  ", @y, "\n";                    # print second sequence
}
close(SEQ2);

# Smith-Waterman algorithm: fill stage

for($i = 0; $i <= $XL; $i++)               # Initialization to zero for the first column
{
  $M[$i][0] = 0;
  $I[$i][0] = 0;
}
for($j = 0; $j <= $YL; $j++)               # Initialization to zero for the first row
{
  $M[0][$j] = 0;
  $I[0][$j] = 0;
}



$endi = 0; $endj = 0;
$Fmax = 0;
$score = 0;

for($i = 1; $i <= $XL; $i++) {                          # Incrementing columns (first sequence)
  for($j = 1; $j <= $YL; $j++) {                        # Recursion, row by row (second sequence)

    $s[$i][$j] = getSubstScore(\@BLOSUM62,$x[$i],$y[$j]);          # Score of a match is computed by the getSubstScore() function (see below)
    $score = $s[$i][$j];


    if ($M[$i-1][$j-1] + $score > $I[$i-1][$j-1] + $score) # fill the M(atch) matrix.
    {
      $M[$i][$j] = $M[$i-1][$j-1] + $score;
    } else {
      $M[$i][$j] = $I[$i-1][$j-1] + $score;
    }

    if ($M[$i][$j] < 0) { $M[$i][$j] = 0; }		# SW: if M(i,j)<0, M(i,j)=0

    if ($M[$i][$j]>$Fmax)
    {
      $Fmax=$M[$i][$j]; $endi=$i; $endj=$j;
    }

    $t1 = $M[$i][$j-1] - $gop;
    $t2 = $M[$i-1][$j] - $gop;
    $t3 = $I[$i][$j-1] - $gep;
    $t4 = $I[$i-1][$j] - $gep;

    $m1 = $t1>$t2?$t1:$t2;
    $m2 = $t3>$t4?$t3:$t4;

    $I[$i][$j] = $m1>$m2?$m1:$m2;

    } # for j
  } # for i

# print filled SW matrix

#print "\nSmith-Waterman score matrix:\n";
#for($i = 1; $i <= $XL; $i++) {
#  print "\n";
#  for($j = 1; $j <= $YL; $j++) {
#    printf "%4i", $M[$i][$j]                            # prints the score matrix (4 characters per number)
#  }
#}

# Smith-Waterman algorithm: traceback stage             # gives the position of the paired bases
                                                        # in each sequence
print "\nTraceback report:\n";
$i = $endi; push @upper, $x[$i];                        # array for sequence 1, from highest score position
$j = $endj; push @lower, $y[$j];                        # array for sequence 2, from highest score position

while($M[$i][$j] > 0)
{
  #printf "\n%4i%4i", $i, $j;                           # prints pair coordinates of alignment path
  # at every cell we have to check whether we arrived here by a match state or an insert state change
  if ($M[$i-1][$j-1] + $s[$i][$j] == $M[$i][$j])        # have a match!!
  {
    push @upper, $x[$i-1];                              # extends alignment array for seq1 (case of match)
    push @lower, $y[$j-1];                              # extends alignment array for seq2
    $i=$i-1; $j=$j-1;                                   # first element of the array is the last of the alignment
  } else {                                              # have an insertion!!
    if($M[$i-1][$j] - $gop == $I[$i][$j] || $I[$i-1][$j] - $gep == $I[$i][$j])   # in seq2!! (either extending gap or starting a new gap)
    {
      push @upper, $x[$i-1];
      push @lower, '-';                                 # gap in seq2
      $i= $i-1;
    } else {                                            # in seq1!!
      push @upper, '-';                                 # gap in seq1
      push @lower, $y[$j-1];
      $j =$j-1;
    }
  }
} # while

# print optimal alignment:

  print "\n\nOptimal alignment score: $Fmax\n";         # maximum score obtained in the matrix

  printf "\nsequence 1: %4i ", $i;                      # $i is start position of alignment for seq1
  while(scalar(@upper) > 0) {print(pop @upper)}         # print seq1 alignment, starting at end of array
  printf "\nsequence 2: %4i ", $j;                      # $j is start position of alignment for seq2
  while(scalar(@lower) > 0) {print(pop @lower)}         # print seq2 alignment, starting at end of array
  print "\n";


exit 0;


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#  this subroutine accepts parameters:
#  1) reference to substitution matrix
#  2) first amino acid
#  3) second amino acid
#  returns:
#  the substitution score between those two residues
#
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
