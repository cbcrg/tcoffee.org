#!/usr/bin/perl
#
# This script reads two fasta-formatted sequence files and performs a basic version of the 
# Smith-Waterman algorithm following approximately the recipe given in Durbin et al. 1998
# "Biological Sequence Analysis - Probabilistic models of proteins and nucleic acids",
# pp. 22-24.
#
$/ = "\n>";                                # entry delimiter for Fasta-formatted sequence library

# input parsing

open(SEQ1, "$ARGV[0]");                    # reading the first sequence
while(<SEQ1>) {                            # major loop over first sequence in the input file
  /^>?([^\n]*)\n(.*)/s;                    # header goes into $1, sequence into $2
  $name = $1;
  $seq = $2;                               # otherwise $2 would be erased during next pattern match
  $seq =~ s/[^\D]//gs;                   # delete all characters in $seq except A, C, G, T.
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
  $seq =~ s/[^\D]//gs;
  @y = split //, " " . $seq;
  $YL = scalar(@y) - 1;                    # $YL is lengt of second sequence
  print "\nSequence 2: $name\n\n";         # print header of second sequence
  print "  ", @y, "\n";                    # print second sequence
  }   

# Smith-Waterman algorithm: fill stage

$r = 1; $q = -2; $d = 0;                  # $r is for a match, $q for a mismatch, $d for a gap 

for($i = 0; $i <= $XL; $i++) {             # Initialization to zero for the first column
  $F[$i][0] = 0
  }
for($j = 0; $j <= $YL; $j++) {             # Initialization to zero for the first row
  $F[0][$j] = 0
  }

$Fmax = 0;  
$endi = 0; $endj = 0;

for($i = 1; $i <= $XL; $i++) {                          # Incrementing columns (first sequence)
  for($j = 1; $j <= $YL; $j++) {                        # Recursion, row by row (second sequence)
    if  ($x[$i] eq $y[$j]) {$s[$i][$j] = $r}            # Score of match is one
    else                   {$s[$i][$j] = $q}            # Score of mismatch is -2
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
  printf "\n%4i%4i", $i, $j;                            # prints pair coordinates of alignment path
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
