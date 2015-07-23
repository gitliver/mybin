#!/bin/bash

echo [start]
echo [pwd] `pwd`
echo [date] `date`

mate1=$1
mate1b=$2
mate2=$3
mate2b=$4
myoutputdir=$5
name=$6
ref=$7

# http://www.linuxjournal.com/content/shell-process-redirection
set +o posix

# bowtie2 -x $ref -1 $mate1 -2 $mate2 -p 2 | samtools view -bS - > ${myoutputdir}/${name}.bam
echo 1
bwa aln -t 4 $ref <( cat $mate1 $mate1b | gunzip ) > ${myoutputdir}/$name.1.sai
echo 2
bwa aln -t 4 $ref <( cat $mate2 $mate2b | gunzip ) > ${myoutputdir}/$name.2.sai
echo 3
bwa sampe $ref ${myoutputdir}/$name.1.sai ${myoutputdir}/$name.2.sai <( cat $mate1 $mate1b | gunzip ) <( cat $mate2 $mate2b | gunzip ) | samtools view -bS - > ${myoutputdir}/${name}.bam

echo [finish]
echo [date] `date`
