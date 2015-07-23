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
ref=/ifs/scratch/c2b2/rr_lab/shares/ref/hg19/faidx/hg19.fa

# map, convert sam to bam, sort bam
# bowtie2 -x $ref -1 $mate1 -2 $mate2 -p 2 | samtools view -bS - | samtools sort - $mylocaloutputdir/${name}.sorted
bowtie2 -x $ref -1 ${mate1},${mate1b} -2 ${mate2},${mate2b} -p 2 | samtools view -bS - > ${myoutputdir}/${name}.bam

echo [finish]
echo [date] `date`

