#!/usr/sbin/perl

#
#   A course on sequence Alignments
# 	Cedric Notredame 2001
#   All rights reserved
#   Files can be redistributed without permission
#   Comercial usage is forbiden
#   Bubble Sort example


# Read the numbers from the command line
@list=@ARGV;

#do the sorting
$ne=$#list;

$swap=1;
while ($swap!=0)
	{
	$swap=0;
	for ($i=0; $i<$ne; $i++)
		{
		if ($list[$i]>$list[$i+1])
			{
			$tmp=$list[$i];
			$list[$i]=$list[$i+1];
			$list[$i+1]=$tmp;
			$swap=1;
			}
		}
	$ne--;
	}
print "@list\n";




                         
