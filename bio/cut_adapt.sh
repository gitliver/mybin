#!/bin/bash

echo [start]
echo [pwd] `pwd`
echo [date] `date`

mate1=$1
mate2=$2
myoutputdir=$3
name=$4

mkdir -p ${myoutputdir}/cutad
myoutputdir=${myoutputdir}/cutad

cutadapt -b GATCGGAAGAGC -B AGATCGGAAGAGCG -o ${myoutputdir}/${name}.R1.cutadapt.fastq -p ${myoutputdir}/${name}.R2.cutadapt.fastq $mate1 $mate2  

gzip ${myoutputdir}/${name}.R1.cutadapt.fastq 
gzip ${myoutputdir}/${name}.R2.cutadapt.fastq 

echo [finish]
echo [date] `date`
