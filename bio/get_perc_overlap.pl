#!/usr/bin/env perl

my $usage = <<_EOUSAGE_;
################################################################################################################
#
# About:
# This script takes a standard tab-delimited blast file and gets the *cumulative* subject and query 
# coverages with bedtools. Because blast is a local aligner, the same query can hit a subject multiple 
# times in different places. It is useful to know the cumulative coverage values and this information 
# cannot be computed by considering rows in isolation.
# 
# The input blast file must be in this format:
# "qseqid sseqid pident alnlength mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
# The sseqid should be in the format, e.g., "gi|164607370|gb|AC193179.2|" because the script presupposes 
# this to pull out the accession, AC193179.2
# 
# The output format is that of Bedtools coverage described here:
# http://bedtools.readthedocs.org/en/latest/content/tools/coverage.html
# I.e., "id feat_start feat_end #_of_features_overlapping feat_overlap_start feat_overlap_end percent_overlap"
#
# Useage:
# Suppose you have a blast file example/blastout.txt. Then: 
# $ cat example/blastout.txt | $0 --prefix mytest --clean
# produces two files: 
# mytest.cov_q.bed
# mytest.cov_s.bed
# where the 7th column of each file is the coverage percent 
#
# Dependencies: 
# Bedtools 
#
# Notes:
# This script goes through a blast file and produces query and subject range bed files, containing all the 
# alignment coordinates in bed format as well as query and suject ref bed files, which are just the full 
# length query and subject coordinates. Then it does a bedtools coverage command to get the coverage.
# To do this, the script makes a tmp bed range file for each unique query-subject pair. It loops thro 
# the blast file and, when the query-subject changes, does the bedtools command and appends the answer
# onto the resultant coverage file.
#
################################################################################################################
_EOUSAGE_

# -----------------------------------------------------------
# PERL REGEX REMINDER

# \d [0-9] Any digit
# \D [^0-9] Any character not a digit
# \w [0-9a-zA-z_] Any "word character"
# \W [^0-9a-zA-z_] Any character not a word character
# \s [ \t\n\r\f] whitespace (space, tab, newline, carriage return, form feed)
# \S [^ \t\n\r\f] Any non-whitespace character

# *      Match 0 or more times
# +      Match 1 or more times
# ?      Match 1 or 0 times
# {n}    Match exactly n times
# {n,}   Match at least n times
# {n,m}  Match at least n but not more than m times
# -----------------------------------------------------------

use strict;
use Data::Dumper;
use Cwd 'abs_path';
use FindBin qw($RealBin);
use Getopt::Long;

my $clean;	
my $prefix;	
my $help;	

GetOptions (	'clean' => \$clean,		# clean up after
		'help' => \$help,		# help 
		'prefix=s' => \$prefix );	# file prefix

if ($help)
{
	print $usage;
	exit;
}

my $bool_first=1;
my $prevqid;    
my $prevsid;  
my $id;
my $counter=1;

my $cmd = "rm -f $prefix.cov_q.bed $prefix.cov_s.bed";
print "$cmd\n";
system($cmd);

while (<STDIN>) 
{

	print $counter,"\n" if ($counter % 1000 == 0 );
	$counter++;

	chomp($_);

	my @line=split("\t",$_); 

	# fields of blast file
	my $qseqid = $line[0];  
	my $sseqid = $line[1];  
	my $pident = $line[2];  
	my $alnlength = $line[3];       
	my $mismatch = $line[4];        
	my $gapopen = $line[5]; 
	my $qstart = $line[6];  
	my $qend = $line[7];    
	my $sstart = $line[8];  
	my $send = $line[9];    
	my $evalue = $line[10]; 
	my $bitscore = $line[11];
	my $qlen = $line[12];
	my $slen = $line[13];

	# split sseqid on "|" to get accession 
	my @a = split(/\|/, $sseqid);
	# accession
	my $sseqid = $a[3];	

	# if first row
	if ($bool_first)
	{
		$prevqid = $qseqid;    
		$prevsid = $sseqid;  
		$bool_first=0;

		$id = $qseqid."@".$sseqid;

		my $str1 = $id."\t"."0"."\t".$slen;
		my $cmd = "echo \"$str1\" > $prefix.sref.bed";
		# print "$cmd\n";
		system($cmd);

		my $str1 = $id."\t"."0"."\t".$qlen;
		my $cmd = "echo \"$str1\" > $prefix.qref.bed";
		# print "$cmd\n";
		system($cmd);

		my $str2;
		if ($sstart <= $send)
		{
			# convert to bedtools/UCSC format, which has the idiotic convention of 0-based start coordinate, end coordinate not included
			$str2 = $id."\t".($sstart-1)."\t".$send;
		}
		else
		{
			$str2 = $id."\t".($send-1)."\t".$sstart;
		}
		my $cmd = "echo \"$str2\" > $prefix.srange.bed";
		# print "$cmd\n";
		system($cmd);
		
		$str2 = $id."\t".($qstart-1)."\t".$qend;
		my $cmd = "echo \"$str2\" > $prefix.qrange.bed";
		# print "$cmd\n";
		system($cmd);
	}
	# if query id and subject id same as previous row
	elsif ($qseqid eq $prevqid && $sseqid eq $prevsid)
	{
		my $str2;
		if ($sstart <= $send)
		{
			$str2 = $id."\t".($sstart-1)."\t".$send;
		}
		else
		{
			$str2 = $id."\t".($send-1)."\t".$sstart;
		}
		my $cmd = "echo \"$str2\" >> $prefix.srange.bed";
		# print "$cmd\n";
		system($cmd);

		$str2 = $id."\t".($qstart-1)."\t".$qend;
		my $cmd = "echo \"$str2\" >> $prefix.qrange.bed";
		# print "$cmd\n";
		system($cmd);
	}
	# if change 
	else
	{
		$id = $qseqid."@".$sseqid;
		$prevqid = $qseqid;    
		$prevsid = $sseqid;  

		$id = $qseqid."@".$sseqid;

		my $str1 = $id."\t"."0"."\t".$slen;
		my $cmd = "echo \"$str1\" >> $prefix.sref.bed";
		# print "$cmd\n";
		system($cmd);

		my $str1 = $id."\t"."0"."\t".$qlen;
		my $cmd = "echo \"$str1\" >> $prefix.qref.bed";
		# print "$cmd\n";
		system($cmd);

		my $str2;
		if ($sstart <= $send)
		{
			$str2 = $id."\t".($sstart-1)."\t".$send;
		}
		else
		{
			$str2 = $id."\t".($send-1)."\t".$sstart;
		}
		my $cmd = "echo \"$str2\" >> $prefix.srange.bed";
		# print "$cmd\n";
		system($cmd);
		
		$str2 = $id."\t".($qstart-1)."\t".$qend;
		my $cmd = "echo \"$str2\" >> $prefix.qrange.bed";
		# print "$cmd\n";
		system($cmd);
	}
}

my $cmd="coverageBed -a $prefix.srange.bed -b $prefix.sref.bed > $prefix.cov_s.bed";
print "$cmd\n";
system($cmd);

my $cmd="coverageBed -a $prefix.qrange.bed -b $prefix.qref.bed > $prefix.cov_q.bed";
print "$cmd\n";
system($cmd);

if ($clean)
{
	system("rm $prefix.srange.bed $prefix.sref.bed $prefix.qrange.bed $prefix.qref.bed");
}
