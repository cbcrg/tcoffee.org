#!/usr/bin/perl
#
# This script reads two fasta-formatted sequence files and performs a basic version of the 
# Smith-Waterman algorithm following approximately the recipe given in Durbin et al. 1998
# "Biological Sequence Analysis - Probabilistic models of proteins and nucleic acids",
# pp. 22-24.
#
$/ = "\n>";                                # entry delimiter for Fasta-formatted sequence library

#
# @x  seq 1
# $XL length of 1st sequence
# @y  seq 2
# $YL length of 2nd sequence
# @s  scoring matrix
# @F  alignment matrix

# input parsing

open(SEQ1, "$ARGV[0]");                    # reading the first sequence
while(<SEQ1>) {                            # major loop over first sequence in the input file
  /^>?([^\n]*)\n(.*)/s;                    # header goes into $1, sequence into $2
  $name = $1;
  $seq = $2;                               # otherwise $2 would be erased during next pattern match
  $seq =~ s/[^ARNDCQEGHILKMFPSTWYV]//gs;   # delete all characters in $seq except amino acids: ARNDCQEGHILKMFPSTWYV
  @x = split //, " " . $seq;               # make sequence array, note $x[0] = ' ', $x[1] = 1st base
  $XL = scalar(@x) - 1;                    # $XL is lengt of first sequence
  print "\nSequence 1: $name\n\n";         # print header of first sequence
  print "  ", @x, "\n";                    # print first sequence
  }

open(SEQ2, "$ARGV[1]");                    # reading the first sequence
while(<SEQ2>) {                            # same as for the first sequence
  /^>?([^\n]*)\n(.*)/s;
  $name = $1;
  $seq = $2;
  $seq =~ s/[^ARNDCQEGHILKMFPSTWYV]//gs;
  @y = split //, " " . $seq;
  $YL = scalar(@y) - 1;                    # $YL is lengt of second sequence
  print "\nSequence 2: $name\n\n";         # print header of second sequence
  print "  ", @y, "\n";                    # print second sequence
  }

# Smith-Waterman algorithm: fill stage

for($i = 0; $i <= $XL; $i++) {             # Initialization to zero for the first column
  $F[$i][0] = 0
  }
for($j = 0; $j <= $YL; $j++) {             # Initialization to zero for the first row
  $F[0][$j] = 0
  }




$gop = -3; $gep = -1;                  # $q for a mismatch, $d for a gap
$d = 0;
$Fmax = 0;
$endi = 0; $endj = 0;

for($i = 1; $i <= $XL; $i++) {                          # Incrementing columns (first sequence)
  for($j = 1; $j <= $YL; $j++) {                        # Recursion, row by row (second sequence)
    $s[$i][$j] = getSubstScore($x[$i],$y[$j]);          # Score of a match is computed by the getSubstScore() function (see below)
    $F[$i][$j] = 0;                                     # zero is the minimum score anyway
    if   ($s[$i][$j] > $F[$i][$j]) {                    # Finding the highest score value
      $F[$i][$j] = $s[$i][$j]                           # for the [i][j] pair. Maximization over the four
      }                                                 # different ways to extend the score matrix,
    if   ($F[$i-1][$j-1] + $s[$i][$j] > $F[$i][$j]) {   # depending on matches, mismatches and gaps.
      $F[$i][$j] = $F[$i-1][$j-1] + $s[$i][$j]
      }
    elsif($F[$i-1][$j]   + $d         > $F[$i][$j]) {
      $F[$i][$j] = $F[$i-1][$j] + $d
      }
    elsif($F[$i][$j-1]   + $d         > $F[$i][$j]) {
      $F[$i][$j] = $F[$i][$j-1] + $d
      }
    if($F[$i][$j] > $Fmax) {                            # change the value of highest score if the new
      $Fmax = $F[$i][$j];                               # one is higher than the old one
      $endi = $i; $endj = $j;                           # keeps track of end positions of highest scoring
      }                                                 # alignment
    }
  }

# print filled SW matrix

print "\nSmith-Waterman score matrix:\n";
for($i = 1; $i <= $XL; $i++) {
  print "\n";
  for($j = 1; $j <= $YL; $j++) {
    printf "%4i", $F[$i][$j]                            # prints the score matrix (4 characters per number)
    }
  }

# Smith-Waterman algorithm: traceback stage             # gives the position of the paired bases
                                                        # in each sequence
print "\nTraceback report:\n";
$i = $endi; push @upper, $x[$i];                        # array for sequence 1, from highest score position
$j = $endj; push @lower, $y[$j];                        # array for sequence 2, from highest score position

while($F[$i][$j] > $s[$i][$j]) {
  #printf "\n%4i%4i", $i, $j;                            # prints pair coordinates of alignment path
  if   ($F[$i-1][$j-1] + $s[$i][$j] == $F[$i][$j]) {
    push @upper, $x[$i-1];                                # extends alignment array for seq1 (case of match)
    push @lower, $y[$j-1];                                # extends alignment array for seq2
    $i=$i-1; $j=$j-1;                                   # first element of the array is the last of the
    }                                                   # alignment
  elsif($F[$i-1][$j]   + $d         == $F[$i][$j]) {
    push @upper, $x[$i-1];
    push @lower, '-';                                   # gap in seq2
    $i= $i-1;
    }
  elsif($F[$i][$j-1]   + $d         == $F[$i][$j]) {
    push @upper, '-';                                   # gap in seq1
    push @lower, $y[$j-1];
    $j =$j-1;
    }
  else {
    die "Error: traceback failed!\n"
    }
  }

# print optimal alignment:

  print "\n\nOptimal alignment score: $Fmax\n";         # maximum score obtained in the matrix

  printf "\nsequence 1: %4i ", $i;                      # $i is start position of alignment for seq1
  while(scalar(@upper) > 0) {print(pop @upper)}         # print seq1 alignment, starting at end of array
  printf "\nsequence 2: %4i ", $j;                      # $j is start position of alignment for seq2
  while(scalar(@lower) > 0) {print(pop @lower)}         # print seq2 alignment, starting at end of array
  print "\n";


exit 0;



sub getSubstScore
{
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




	$aaOrder = "ARNDCQEGHILKMFPSTWYV";
	$aa1 = shift;
	$aa2 = shift;

	$pos1 = index($aaOrder, $aa1);
	$pos2 = index($aaOrder, $aa2);

	return ( $BLOSUM62[$pos1][$pos2] );
}
