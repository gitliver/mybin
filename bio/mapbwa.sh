#!/bin/bash

echo [start]
echo [pwd] `pwd`
echo [date] `date`

mate1=$1
mate2=$2
myoutputdir=$3
name=$4
ref=$5

# http://www.linuxjournal.com/content/shell-process-redirection
set +o posix

# bowtie2 -x $ref -1 $mate1 -2 $mate2 -p 2 | samtools view -bS - > ${myoutputdir}/${name}.bam
echo 1
bwa aln -t 4 $ref <( gunzip --stdout $mate1 ) > ${myoutputdir}/$name.1.sai
echo 2
bwa aln -t 4 $ref <( gunzip --stdout $mate2 ) > ${myoutputdir}/$name.2.sai
echo 3
bwa sampe $ref ${myoutputdir}/$name.1.sai ${myoutputdir}/$name.2.sai <( gunzip --stdout $mate1 ) <( gunzip --stdout $mate2 ) | samtools view -bS - > ${myoutputdir}/${name}.bam

echo [finish]
echo [date] `date`

