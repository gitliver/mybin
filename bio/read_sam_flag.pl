#!/usr/bin/env perl

# About:
# script to interpret the samtools flag

# Useage:
# Suppose you want to see what bwa flag 99 is:
# echo 99 | read_sam_flag.pl 

# Notes:
# Refer to http://samtools.sourceforge.net/samtools.shtml for more info about the flag

while (<STDIN>)
{
	chomp $_; 
	
	$binarynum=sprintf ("%b", $_);	

#	print $binarynum,"\n";
#	print $_,"=",$binarynum,"\n";
	
	$rbinarynum=reverse($binarynum);

#	print length($rbinarynum),"\n";

	$padding="";

	for ($count = 1; $count <= 11-length($rbinarynum); $count++) 
	{
 		$padding=$padding."0";
 	}
 	
# 	print $padding,"\n";

	$rbinarynum=$rbinarynum.$padding;
		
	print substr($rbinarynum, 0, 1),"\t","the read is paired in sequencing","\n";
	print substr($rbinarynum, 1, 1),"\t","the read is mapped in a proper pair","\n";
	print substr($rbinarynum, 2, 1),"\t","the query sequence itself is unmapped","\n";
	print substr($rbinarynum, 3, 1),"\t","the mate is unmapped","\n";
	print substr($rbinarynum, 4, 1),"\t","strand of the query (1 for reverse)","\n";
	print substr($rbinarynum, 5, 1),"\t","strand of the mate","\n";
	print substr($rbinarynum, 6, 1),"\t","the read is the first read in a pair","\n";
	print substr($rbinarynum, 7, 1),"\t","the read is the second read in a pair","\n";
	print substr($rbinarynum, 8, 1),"\t","the alignment is not primary","\n";
	print substr($rbinarynum, 9, 1),"\t","QC failure","\n";
	print substr($rbinarynum, 10, 1),"\t","optical or PCR duplicate","\n";
	
	print "\n";		
	
}
