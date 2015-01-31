#!/usr/bin/env perl

# merge file on rows, using the first field of the first column as a key

use Data::Dumper; 
use strict; 

my %h=();

while (<STDIN>)
{
	chomp($_); 
	my @a=split; 

	if ( not $h{$a[0]} ) 
	{ 
		$h{$a[0]} = $a[1]; 
	} 
	else 
	{ 
		$h{$a[0]} = $h{$a[0]}.",".$a[1]; 
	}
}

foreach my $key ( keys %h ) 
{
	print $key,"\t",$h{$key},"\n"
}

# print Dumper \%h
