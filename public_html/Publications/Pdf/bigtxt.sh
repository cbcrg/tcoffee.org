#!/bin/bash

RESULT="all.txt"

if [ -e $RESULT ]; then 
  rm -rf $RESULT
fi 

for FILE in `ls *.pdf`
do 
  echo Processing $FILE ... 
  pdftotext $FILE /dev/stdout | cat >> $RESULT
done
