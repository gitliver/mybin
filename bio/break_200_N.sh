#!/bin/bash

myhelp=$(cat <<_EOF_

Take a fasta file and break it into smaller entries on stretches of 200, or more, consequtive "N" characters. 
Also, all Ns at beginning and end of a line are replaced.
Entries in the output fasta are given uniq numeric IDs.

Usage:
break_200_N.sh <file>

Dependencies:

fastajoinlines

_EOF_)

if [  $# == 0 -o "$1" == "-h" -o "$1" == "-help" -o "$1" == "--help" ]; then
	shift; echo "$myhelp"; exit;
fi

filein=$1
cutoff=5 # print if resultant sequence > cutoff

# sloppy - shouldn't really combine awk and perl!

cat $filein | fastajoinlines | awk '$0 !~ />/' | perl -pe 's/^N+//;s/N+$//;s/N{200,}/\n/g;' | awk -v cutoff=$cutoff 'BEGIN{counter=1}{if (length >= cutoff) {print ">"counter; print; counter++}}'
