#!/bin/bash

echo [start]
echo [pwd] `pwd`
echo [date] `date`

mate1=$1
mate2=$2
myoutputdir=$3
name=$4
ref=$5

# map, convert sam to bam, sort bam
# bowtie2 -x $ref -1 $mate1 -2 $mate2 -p 2 | samtools view -bS - | samtools sort - $mylocaloutputdir/${name}.sorted
bowtie2 -x $ref -1 $mate1 -2 $mate2 -p 2 | samtools view -bS - > ${myoutputdir}/${name}.bam

echo [finish]
echo [date] `date`
