#!/usr/bin/env perl

use strict;
use Data::Dumper;

my $usage = <<_EOUSAGE_;

################################################################################################################
#
#  About:
#
#  Input a blast file with degenerate qseqid's and this scirpt 
#  returns a blast file with only top hits --- i.e., only uniq qseqid's
#
#  Useage example:
#  $0 input.blast 
#
#  example:
#
#  if the input were:
#
#  4	gi|292384466|gb|GU938867.1|	99.34	152	1	0	1	152	7	158	4e-71	 276	152	591
#  4	gi|265524965|gb|GU071091.1|	99.34	152	1	0	1	152	8173	8324	4e-71	 276	152	39778
#  4	gi|37956638|gb|AY264774.1|	99.34	152	1	0	1	152	8173	8324	4e-71	 276	152	39938
#  6	gi|265524965|gb|GU071091.1|	100.00	150	0	0	1	150	2886	2737	2e-71	 278	300	39778
#  6	gi|265524965|gb|GU071091.1|	98.08	52	1	0	140	191	578	527	3e-15	91.6	300	39778
#
#  the output would be:
#
#  4	gi|292384466|gb|GU938867.1|	99.34	152	1	0	1	152	7	158	4e-71	 276	152	591
#  6	gi|265524965|gb|GU071091.1|	100.00	150	0	0	1	150	2886	2737	2e-71	 278	300	39778
####################################################################################################################

_EOUSAGE_

#my $arg1 = shift;
#my $arg2 = shift;

my $arg1 = $ARGV[0];

if ($arg1 eq "-h" or $arg1 eq "--h" or $arg1 eq "-help" or $arg1 eq "--help" or scalar(@ARGV) == 0)
{
	print $usage;
	exit;
}

open my $infile, '<', $arg1;
#open my $outfile, '>', $arg2;				

# hash of qid's
my %tmph=();

# ASSUME qid is first column

while (<$infile>)
{	
	my @ln=split; 
	
	# ASSUME best hit is first
	# if no preexisting entry in the hash (i.e., if new qid), print line 
	if (!($tmph{$ln[0]})) 
	{
		#print $outfile $_;
		print $_;
	}
	
	# hash qid to 1  
	$tmph{$ln[0]}=1;			
}
						 
close($infile); 
