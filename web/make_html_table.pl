#!/usr/bin/env perl

# About:
# Make an HTML table from stuff piped in from stdin 

# Example:
# Suppose
# $ cat example/file_table.txt
# 1     a
# 2     b
# 3     c

# Useage: 
# $ cat example/file_table.txt | ./make_html_table.pl 
# <!DOCTYPE html>
# <html>
# <table border="1">
# <tr><td>1</td><td>a</td></tr>
# <tr><td>2</td><td>b</td></tr>
# <tr><td>3</td><td>c</td></tr>
# </table>
# </html>

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
use warnings;
use Data::Dumper;
use Getopt::Long;
use FindBin qw($RealBin);
use Cwd 'abs_path';

my $usage = <<_EOUSAGE_;

################################################################################################################
#
#  Produce a .html file with a table from tab-delimited data 
#  
#  Useage example:
#  $0 --input data.txt --title "My Table"
#  cat data.txt | $0 
#
#  Required inputs:
#
#  Optional inputs:
#
#  --input <string>			input file 
#
#  --title <string>			Title in .html file
#
#  --outputdir <string>			Output directory (default: cwd)
#
#  --help				help
#
####################################################################################################################

_EOUSAGE_

### VARIABLES
my $infile;						# input file
my $title;						# title
my $help;						# help bool
my $outputdir=Cwd::getcwd();				# output dir (default: cwd)
my $numarg=scalar(@ARGV);				# number of args

### OPTIONS
GetOptions (	'input=s' => \$infile,			# input file
		'title=s' => \$title,	        	# title
		'outputdir=s' => \$outputdir,	        # outputdir
		'help' => \$help);		        # bool_help

# if ( $help || $numarg == 1 || (not defined($infile) ) ) {print $usage; exit;}
if ( $help ) {print $usage; exit;}

### MAIN
print "<!DOCTYPE html>\n";
print "<html>\n";

if ($title)
{
	print "<h2>$title</h2>\n";
}


print "<table border=\"1\">\n";
while (<>)
{
        chomp($_); 
        my @line=split("\t",$_); 

        print "<tr>";
        foreach my $elt (@line)
        {
                print "<td>$elt</td>";
        }
        print "</tr>\n";
}
print "</table>\n";
print "</html>\n";
