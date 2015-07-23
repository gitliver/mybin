#!/bin/bash

echo [start]
echo [pwd] `pwd`
echo [date] `date`

bam=$1
name=$2
outdir=$3

cd $outdir

# sort bam
samtools sort ${bam} ${name}.sorted
samtools index ${name}.sorted.bam

echo [finish]
echo [date] `date`
