#!/bin/bash

# About:
# format lyrics (copied mostly from http://www.azlyrics.com) for LaTex, 
# with the goal of making a songbook
# Specifically, add LaTeX style carriage returns where there are 
# real carriage returns and linebreaks where there are spaces

# Example:
# Suppose that 
# $ cat old_apartment_lyrics.txt 
# Broke into the old apartment
# This is where we used to live
# Broken glass, broke and hungry
# Broken hearts and broken bones
# This is where we used to live

# Useage:
# $ ./format_lyrics_4_tex.sh --name "The Old Apartment" --author "Barenaked Ladies" --input example/old_apartment_lyrics.txt 

# The outcome is:
# \section*{Barenaked Ladies}
# \addcontentsline{toc}{section}{Barenaked Ladies}
# 
# \subsection*{The Old Apartment}
# \addcontentsline{toc}{subsection}{The Old Apartment}
# 
# Broke into the old apartment\\
# This is where we used to live\\
# Broken glass, broke and hungry\\
# Broken hearts and broken bones\\
# This is where we used to live
# 
# \newpage

# help section
myhelp=$( cat <<_EOF_

format lyrics for latex, by preserving carriage returns as they appear in the text 

Usage:
./format_lyrics_4_tex.sh [ --name <string> ] [ --author <string> ] [ --input <file> ]

Example: Suppose

# $ cat \$input 
# hello
# world
#
# hello

Then running this script will produce:

# hello\\\\ 
# world\\linebreak   
#
# hello\\newpage 

_EOF_ )

# if no arguments, print help and exit
if [ $# == 0 ]; then echo "$myhelp"; exit; fi

# get options 
while [ $# -gt 0 ]; do

        if [  "$1" == "-h" -o "$1" == "-help" -o "$1" == "--help" ]; then
                shift; echo "$myhelp"; exit; 

        elif [  "$1" == "-name" -o "$1" == "--name" ]; then
                shift; myname=$1; shift

        elif [  "$1" == "-author" -o "$1" == "--author" ]; then
                shift; myauthor=$1; shift

        elif [  "$1" == "-input" -o "$1" == "--input" ]; then
                shift; myinput=$1; shift

        fi

done

# print LaTeX section and subsection  
if [ "$myauthor" ]; then
	echo "\section*{$myauthor}"
	echo "\addcontentsline{toc}{section}{$myauthor}"
	echo
fi

if [ "$myname" ]; then
	echo "\subsection*{$myname}"
	echo "\addcontentsline{toc}{subsection}{$myname}"
	echo
fi

# if file exists, add LaTeX newlines and linebreaks to preserve current format
if [ -e "$myinput" ]; then

	cat $myinput | awk '
	{
		# if first line, store it in variable "previous"
		if (NR==1)
		{
			prev=$0
		}
		else
		{
			# if line empty, add a linebreak to previous line
			if ($0 == "")
			{
				print prev"\\linebreak"; 
			}
			else if (prev == "") # else if prev empty just print it 
			{
				print prev;
			}
			else # otherwise, add an ordinary carrige return
			{
				print prev"\\\\";
			}
			
			prev=$0
		}	
	} 
	END{
		# we are one line behind, so print last line in end block
		# we hope prev is not empty
		if (prev == "") 
		{
			print prev;
		}
		else # for last one, add newpage
		{
			print prev;
			printf "\n";
			print "\\newpage";
			printf "\n";
		}
	}' 

fi
