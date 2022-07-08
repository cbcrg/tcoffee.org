#!/usr/bin/env perl
use Env;

$cl=join (" ",@ARGV);

if (($cl=~/\-dir=(\S+)/))
  {
    $dir=$1;
  }
else
  {
    $dir="./";
  }

if (($cl=~/\-age=(\S+)/))
  {
    $max_age=$1;
  }
else
  {
    $max_age=190;
  }

if ($cl=~/\-size=(\S+)/)
  {
    $max_size=$1;
  }
else
  {
    $max_size=100;
  }

$max_size*=1000000;

if ( ! -d $dir)
  {
    exit (EXIT_FAILURE);
  }


@list=`ls -t1 $dir`;
foreach $f (@list)
  {
    chomp $f;
    $f="$dir/$f";
    
    $s= -s ($f);
    $a= -M ($f);
    $tot_size+=$s;
    $TOT+=$s;
    
    if ($tot_size >=$max_size || $a>=$max_age)
      {
	unlink ($f);
	$remove++;
	$tot_size-=$s;
      }
    else
      {
	$keep++;
      }
  }

exit (EXIT_SUCCESS);



