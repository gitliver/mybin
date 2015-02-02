#!/usr/bin/env perl

use strict;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use Cwd 'abs_path';
use lib "~/SCRIPTS/perl/BioPerl-1.6.1/Bio";
use Bio::Seq;

use FindBin qw($RealBin);

my $usage = <<_EOUSAGE_;

################################################################################################################
#
#  Produce kmers from an amino acid sequence; translate kmers into nucleotides if specified
#  
#  Useage example:
#  cat my_alignment.fasta | $0 --name sample --kmer 50 --slide 10 --codontable aa2nt.txt --header
#
#  Required inputs:
#  --name <string>				name 
#  --kmer <int>					aa kmer length  
#  --slide <int>				aa slide length 
#
#  Optional inputs:
#  --codontable					tab-delimited file of col1=aa, col2=codon (*** rna codon - use "U" not "T" ***)  
#  --header					print header
#  --verbose					verbose mode
#
#  Dependencies: 
#  BioPerl
#
#  Example: 
#  cat VP40-aa.fasta | $0 --name VP40-aa --kmer 100 --slide 10
#
#  Example output: 
#  name	start_postion	#uniq_seqs_at_start_position	aa-mer
#  VP40-aa	0	2	MFSKHVTLPPPPYNPSSPEGLYYNPDAGKKHGNPHSHPVPHTEHSRASNTIRPTADFSLDYDSSSASGAISAFMLEAYVNVISNNKVLLKLVPLWLPLGV
#  VP40-aa	0	2	MRRGVLPTAPPAYNDIAYSMSILPTRPSVIVNETKSDVLAVPGADVPSNSMRPVADDNIDHSSHTPSGVASAFILEATVNVISGTKVLMKQIPIWLPLGV
#  VP40-aa	10	2	PAYNDIAYSMSILPTRPSVIVNETKSDVLAVPGADVPSNSMRPVADDNIDHSSHTPSGVASAFILEATVNVISGTKVLMKQIPIWLPLGVADQKIYSFDS
#  VP40-aa	10	2	PPYNPSSPEGLYYNPDAGKKHGNPHSHPVPHTEHSRASNTIRPTADFSLDYDSSSASGAISAFMLEAYVNVISNNKVLLKLVPLWLPLGVAGQDLYSFDS
#  VP40-aa	20	2	LYYNPDAGKKHGNPHSHPVPHTEHSRASNTIRPTADFSLDYDSSSASGAISAFMLEAYVNVISNNKVLLKLVPLWLPLGVAGQDLYSFDSTASALLIASY
#
####################################################################################################################

_EOUSAGE_

my $infile;
my $kmer;
my $slide;
my $codontable;
my $header;
my $verbose;

GetOptions (	'name=s' => \$infile,			# file name
		'kmer=s' => \$kmer,			# aa kmer length
		'slide=s' => \$slide,			# slide length (in aa base pairs)
		'codontable=s' => \$codontable,		# codon table (if this one is included, the aa-mers will be translated to nt-mers)
		'verbose' => \$verbose,			# verbose bool
		'header' => \$header);			# header bool - if on, print header
            
die print $usage if (!( defined($infile) && defined($kmer) && defined($slide) ));

my %h = ();		# two-tiered hash of global alignment start position and fasta ID to kmer	
my ($id);		# fasta ID

# e.g., with slide==1, kmer==5 and this fasta file:
# >ebolavirus
# MFSKHVTLPPP
# >othervirus
# YFSKHVTLPPP
# you get:
# h{0}{ebolavirus} --> MFSKH
# h{0}{othervirus} --> YFSKH
# h{1}{ebolavirus} --> FSKHV
# etc

my %hend = ();		# two-tiered hash of global alignment start position and fasta ID to global alignment end position
			# BAD STYLE ! the proper technique is a three-tiered hash sample start end -> kmer

my %htot = ();		# a three-tiered hash sample start end -> kmer

# make hash
while (<STDIN>)
{
	chomp $_;
	
	if ( $_ =~ m/^>/)
	{
		$id = substr($_,1);
	}
	elsif ( not m/^(\s*)$/ ) 
	{		
		my $offset=0;
                                      
		# line is the whole line
		my $seq = $_;
		# print $seq,"\n";
		my $seqlen=length($seq);
		
		# if offset < seq length (0 based counting), loop thro string
		while ($offset < $seqlen) 
		{
			my $res="";
			
			# OLD WAY:
			# get substring of length $kmer starting from position $offset (0 based counting)
			if (0)
			{
				$res=substr($seq,$offset,$kmer); 
			}
			
			# the above was the OLD way of getting the kmer. However, the new directive is, get the kmer jumping over indels:
			
			# NEW WAY:
			# if start on indel, skip
			if ( substr($seq, $offset, 1) ne "-" )
			{
				
				my $counter=0;  # counter counts non indel (dash) bases
				my $counter2=0; # counter2 moves over the sequence
				
				# loop until counter reaches kmer OR the end of the sequence is hit
				while ( $counter < $kmer && ($offset + $counter2) < $seqlen ) 
				{
					# print "counter $counter \n";
					# print "o+c2 ",($offset + $counter2)," \n";
					
					# get base 
					my $base=substr($seq, $offset + $counter2, 1);
	
					# if base not equals dash
					if ( $base ne "-" )
					{	
						$res=$res.$base;		 
						$counter++;
					}
					
					$counter2++;
	
				}
				
				print "id: $id, kmer: ",$res,"\n" if ($verbose);
				
				# only print if the kmer length is as specified (we dont want any substrings smaller than kmer)
				if (length($res)==$kmer)
				{
					# hash {start position} {id} --> kmer
					$h{$offset}{$id}=$res;

					# hash {start position} {end position} {id} --> kmer
					$htot{$offset}{($offset + $counter2 - 1)}{$id}=$res;
				}
			
			}
			
			# increase the offset by $slide
			$offset = $offset + $slide;
		}
	}       
}

print "\nhash(start)(end)(id) -> kmer \n" if ($verbose);
print Dumper \ %htot if ($verbose);

if (0)
{
# loop through outer keys
foreach my $key ( sort keys %h )
{
	print "key: $key \n";

	# get the inner hash
	my %h2 = %{$h{$key}};
	
	foreach my $key2 ( keys %h2 )
	{
		print "key2: $key2, value: $h2{$key2} \n";		
	}
}
}

# old block using %h ... redo below using %htot (i.e., also including end position)
if (0)
{
# print header if option
my $filehead="filename"."\t"."start_postion"."\t"."#seqs"."\t"."#uniq_seqs_at_start_position"."\t"."aa-mer"."\n";
$filehead="filename"."\t"."start_postion"."\t"."#seqs"."\t"."#uniq_seqs_at_start_position"."\t"."aa-mer"."\t"."nt-mer"."\n" if ($codontable);
print $filehead if ($header);

# loop through outer keys (the alignment postion)
foreach my $key ( sort keys %h )
{
	# print "key: $key \n";

	# get the inner hash (the fasta file entry)
	my %h2 = %{$h{$key}};
	
	# print $key,"\t";
	# print Dumper \ %h2;

	# make a hash to hold uniq kmers
	my %h3 = ();

	# make var to count uniq kmers
	my $number_uniq_kmer = 0;
	
	# loop thro inner hash (fasta entries) and add UNIQUE kmers as key entries to %h3 (hashed to number of times the kmer appears)
	foreach my $key2 ( keys %h2 )
	{
		
		# print "key2: $key2, val: $h2{$key2}, num: $number_uniq_kmer \n";		
		
		if ( ! $h3{$h2{$key2}} )
		{			
			# dont include stop codons, indels
			if ( $h2{$key2} !~ m/\*/ && $h2{$key2} !~ m/\-/ && $h2{$key2} !~ m/\?/ )
			{
				$h3{$h2{$key2}}=1;
				$number_uniq_kmer++;
			}
		}
		else
		{
			# count how many times kmer occurs
			$h3{$h2{$key2}}++;
		}		
	}

	# loop thro UNIQUE kmers
	foreach my $key3 ( keys %h3 )
	{
		# if codons provided, translate AA 2 NT
		if ($codontable)
		{
			my $href = col1_to_col2_hash($codontable);
			my %codonhash=%$href;
			
			my $ntseq="";	# NT sequence
			# loop thro AA and convert AA 2 NT
			for (my $i = 0; $i < length($key3); $i++) 
			{
				# get AA
				my $letter=substr($key3,$i,1);
				
				# if entry in hash, translate
				if ($codonhash{$letter})
				{
					$ntseq=$ntseq.$codonhash{$letter};
				}
				else
				{
					die "[error] codon translation not found for $letter \n";
				}	
			}
			
			# convert RNA -> DNA (i.e., change U to T) and print fasta in a file sans leftstring and rightstring. This will serve as the bowtie2 reference
			$ntseq =~ s/U/T/g;			
			
			# print inputfile, globalalignposition, #kmer, diversity, seq, ntseq	
			print $infile,"\t",$key,"\t",$h3{$key3},"\t",$number_uniq_kmer,"\t",$key3,"\t",$ntseq,"\n";
		}
		else
		{
			print $infile,"\t",$key,"\t",$h3{$key3},"\t",$number_uniq_kmer,"\t",$key3,"\n";
		}
	}
	
}
}

# print header if option
my $filehead="filename"."\t"."start_pos"."\t"."end_pos"."\t"."#seqs_at_start_end"."\t"."#uniq_seqs_at_start_end"."\t"."aa-mer"."\n";
$filehead="filename"."\t"."start_pos"."\t"."end_pos"."\t"."#seqs_at_start_end"."\t"."#uniq_seqs_at_start_end"."\t"."aa-mer"."\t"."nt-mer"."\n" if ($codontable);
print $filehead if ($header);

# loop through outer keys - the alignment start postion
foreach my $key ( sort keys %htot )
{
	# get the inner hash (the fasta file entry)
	my %h2 = %{$htot{$key}};
	
	print "pos_start: $key\n" if ($verbose);
	# print Dumper \ %h2;	
	
	# loop through outer keys - the alignment end postion
	foreach my $key2 ( sort keys %h2 )
	{
		# get the inner hash (the fasta file entry)
		my %h3 = %{$h2{$key2}};
		
		print "pos_end: $key2\n" if ($verbose);
		print Dumper \ %h3 if ($verbose);
		
		# make a hash to hold uniq kmers
		my %h_kmer = ();
	
		# make var to count uniq kmers
		my $number_uniq_kmer = 0;
		
		# loop thro inner hash (fasta entries) and add UNIQUE kmers as key entries to %h_kmer (hashed to number of times the kmer appears)
		foreach my $key3 ( keys %h3 )
		{
			# print "key2: $key2, val: $h2{$key2}, num: $number_uniq_kmer \n";		
			
			# get the kmer sequence 
			my $kmer_seq = $h3{$key3};
			
			if ( ! $h_kmer{$kmer_seq} )
			{			
				# dont include stop codons, indels
				if ( $kmer_seq !~ m/\*/ && $kmer_seq !~ m/\-/ && $kmer_seq !~ m/\?/ )
				{
					$h_kmer{$kmer_seq}=1;
					$number_uniq_kmer++;
				}
			}
			else
			{
				# count how many times kmer occurs at a given start and end
				$h_kmer{$kmer_seq}++;
			}		
		}

		# loop thro UNIQUE kmers for given start and end position
		foreach my $key4 ( keys %h_kmer )
		{
			# if codons provided, translate AA 2 NT
			if ($codontable)
			{
				my $href = col1_to_col2_hash($codontable);
				my %codonhash=%$href;
				
				my $ntseq="";	# NT sequence
				# loop thro AA and convert AA 2 NT
				for (my $i = 0; $i < length($key4); $i++) 
				{
					# get AA
					my $letter=substr($key4,$i,1);
					
					# if entry in hash, translate
					if ($codonhash{$letter})
					{
						$ntseq=$ntseq.$codonhash{$letter};
					}
					else
					{
						die "[error] codon translation not found for $letter \n";
					}	
				}
				
				# convert RNA -> DNA (i.e., change U to T) and print fasta in a file sans leftstring and rightstring. This will serve as the bowtie2 reference
				$ntseq =~ s/U/T/g;			
				
				# print inputfile, align_start_position, align_end_position, #kmer at given start and end, diversity, aaseq, ntseq	
				print $infile,"\t",$key,"\t",$key2,"\t",$h_kmer{$key4},"\t",$number_uniq_kmer,"\t",$key4,"\t",$ntseq,"\n" if not ($verbose);
			}
			else
			{
				# print inputfile, align_start_position, align_end_position, #kmer at given start and end, diversity, aaseq	
				print $infile,"\t",$key,"\t",$key2,"\t",$h_kmer{$key4},"\t",$number_uniq_kmer,"\t",$key4,"\n" if not ($verbose);
			}
		}				
	}
}

# return hash of column1 to column2 (easily done with the map function)
sub col1_to_col2_hash  
{		
	# When the subroutine is called any parameters are passed as a list in the special @_ (array form of $_)
	# @_ exists magically - anything passed in

	# arg 1 - file name

	my $infile = shift; 
	
	my %h = ();				# hash	

	if ( -s $infile )		# if file nonzero
	{	
		open(my $fh, '<', $infile);
		
		my $skipline = 0; 		# skip first line (header) ... this is false now, so the first line wont be skipped
			
		# loop thro file
		while (<$fh>)
		{
			if (!$skipline)
			{
				chomp $_;
					
				# line is the whole line
				my $line = $_;
				
				# split on ws
				my @row = split;
				
				$h{$row[0]} = $row[1];
			}	
			$skipline = 0;
		}
		
			close($fh);
	}

	return \%h;
}
