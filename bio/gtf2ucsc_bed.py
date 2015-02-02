#!/usr/bin/env python

import sys
import re

# Learning Python! 

# A script which converts exon ranges in gtf format into UCSC appropriate bed format (as specified here: https://genome.ucsc.edu/FAQ/FAQformat.html#format1)
# Useage: 
# gtf2ucsc_bed.py my.gtf
# Verbose Mode:
# gtf2ucsc_bed.py my.gtf 1

# example input:

# 5     USAMRIID        exon    1629    3068    .       -       .       gene_id "XLOC_033544"; transcript_id "TCONS_00092538"; exon_number "1"; gene_name "POLR3K"
# 5     USAMRIID        exon    7918    8005    .       -       .       gene_id "XLOC_033544"; transcript_id "TCONS_00092538"; exon_number "2"; gene_name "POLR3K"
# 5     USAMRIID        exon    10109   10276   .       -       .       gene_id "XLOC_033544"; transcript_id "TCONS_00092538"; exon_number "3"; gene_name "POLR3K"
# 5     USAMRIID        exon    7268    7291    .       -       .       gene_id "XLOC_033544"; transcript_id "TCONS_00092539"; exon_number "1"; gene_name "POLR3K"
# 5     USAMRIID        exon    7918    8005    .       -       .       gene_id "XLOC_033544"; transcript_id "TCONS_00092539"; exon_number "2"; gene_name "POLR3K"
# 5     USAMRIID        exon    10109   12778   .       -       .       gene_id "XLOC_033544"; transcript_id "TCONS_00092539"; exon_number "3"; gene_name "POLR3K"
# 10    USAMRIID        exon    2669    4520    .       +       .       gene_id "XLOC_002358"; transcript_id "TCONS_00006620"; exon_number "1"; gene_name "LOC721591"
# 10    USAMRIID        exon    4981    5890    .       +       .       gene_id "XLOC_002358"; transcript_id "TCONS_00006620"; exon_number "2"; gene_name "LOC721591"
# 10    USAMRIID        exon    2669    4511    .       +       .       gene_id "XLOC_002358"; transcript_id "TCONS_00006621"; exon_number "1"; gene_name "LOC721591"
# 10    USAMRIID        exon    5149    5890    .       +       .       gene_id "XLOC_002358"; transcript_id "TCONS_00006621"; exon_number "2"; gene_name "LOC721591"
# 10    USAMRIID        exon    2678    4148    .       +       .       gene_id "XLOC_002358"; transcript_id "TCONS_00006624"; exon_number "1"; gene_name "LOC721591"
# 10    USAMRIID        exon    4594    4678    .       +       .       gene_id "XLOC_002358"; transcript_id "TCONS_00006624"; exon_number "2"; gene_name "LOC721591"
# 13    USAMRIID        exon    3576    7683    .       -       .       gene_id "XLOC_008166"; transcript_id "TCONS_00023412"; exon_number "1"; gene_name "FILIP1"
# 13    USAMRIID        exon    21180   21237   .       -       .       gene_id "XLOC_008166"; transcript_id "TCONS_00023412"; exon_number "2"; gene_name "FILIP1"
# 13    USAMRIID        exon    24379   27184   .       -       .       gene_id "XLOC_008166"; transcript_id "TCONS_00023412"; exon_number "3"; gene_name "FILIP1"
# 13    USAMRIID        exon    67096   67274   .       -       .       gene_id "XLOC_008166"; transcript_id "TCONS_00023412"; exon_number "4"; gene_name "FILIP1"
# 13    USAMRIID        exon    76511   76684   .       -       .       gene_id "XLOC_008166"; transcript_id "TCONS_00023412"; exon_number "5"; gene_name "FILIP1"
# 13    USAMRIID        exon    127251  127532  .       -       .       gene_id "XLOC_008166"; transcript_id "TCONS_00023412"; exon_number "6"; gene_name "FILIP1"
# 13    USAMRIID        exon    203289  203742  .       -       .       gene_id "XLOC_008166"; transcript_id "TCONS_00023412"; exon_number "7"; gene_name "FILIP1"
# 13    USAMRIID        exon    3576    7683    .       -       .       gene_id "XLOC_008166"; transcript_id "TCONS_00023413"; exon_number "1"; gene_name "FILIP1"

# example output:

# 10    2668    5890    LOC721591.1     1000    +       .       .       0       2       1852,910        2668,4980
# 10    2668    5890    LOC721591.2     1000    +       .       .       0       2       1843,742        2668,5148
# 10    2677    4678    LOC721591.3     1000    +       .       .       0       2       1471,85 2677,4593
# 13    3575    203742  FILIP1.1        1000    -       .       .       0       7       4108,58,2806,179,174,282,454    3575,21179,24378,67095,76510,127250,203288
# 13    3575    7683    FILIP1.2        1000    -       .       .       0       1       4108    3575
# 5     7267    12778   POLR3K.1        1000    -       .       .       0       3       24,88,2670      7267,7917,10108
# 5     1628    10276   POLR3K.2        1000    -       .       .       0       3       1440,88,168     1628,7917,10108

# define dicts
dtran={}  # transcript id to chr, gene, strand mapping
dgene={}  # gene to transcript mapping
dstlen={} # transcript id to start position, and length mapping

# verbose
verbose = 0

if (len(sys.argv) == 1):
    print('''A script which converts exon ranges in gtf format into UCSC appropriate bed format (as specified here: https://genome.ucsc.edu/FAQ/FAQformat.html#format1).
The ninth column is assumed to be of the format gene_id "X"; transcript_id "X"; exon_number "X"; gene_name "X"

Useage: 
gtf2ucsc_bed.py my.gtf

Verbose Mode:
gtf2ucsc_bed.py my.gtf 1
''');
    exit(0)
if (len(sys.argv) > 2):
    verbose = int(sys.argv[2])

with open(sys.argv[1], "r") as f:
    for line in f:
	# break line on tabs
	listcols = line.split("\t")

	# get columns
	chr = listcols[0]
	exstart = int(listcols[3])
	exend = int(listcols[4])
	strand = listcols[6]
	attr = listcols[8]
	# print(attr),

	# extract strings from attributes column 
	match = re.search(r'gene_id "(\S+)"; transcript_id "(\S+)"; exon_number "(\S+)"; gene_name "(\S+)"', attr)

	if match:
	    # get name
	    geneid = match.group(1)
	    tranid = match.group(2)
	    exonnum = match.group(3)
	    genename = match.group(4)

	    # map transcript isoform id to gene
	    dtran[tranid] = (chr,genename,strand)

	    # map gene to transcript isoform ids
	    if (genename in dgene):
	        dgene[genename][tranid] = 1
	    else:
	        dgene[genename] = {}

	    # map transcript isoform id to its exon starts and lengths
	    if (tranid in dstlen):
	        dstlen[tranid].append((exstart, exend - exstart + 1))
	    else:
	        dstlen[tranid] = [(exstart, exend - exstart + 1)]


if verbose:
    print(dtran)
    print
    print(dgene)
    print
    print(dstlen)
    print

# keys are genes
for mygene in dgene:
	
    if verbose:
        print("***" + mygene + "***")

    # define isoform num to iterate over
    isoform_number = 1 

    for myisoform in dgene[mygene]:
        newgenename=mygene + "." + str(isoform_number)
        isoform_number += 1

	startlist = ""
	endlist = ""
	lenlist = ""
	start_coord = 0
	end_coord = 0
	num_blocks = len(dstlen[myisoform])

	# boolean to track if we're at the first exon in the isoform
	isfirst = 1

	# get starts and lengths 
        for mytup in dstlen[myisoform]:
		# convert format to idiotic UCSC convention: must subtract one from the start
	    if isfirst:
                start_coord = mytup[0] - 1
	        isfirst = 0
	    startlist += "," + str(mytup[0] - 1 - start_coord) # this is an OFFSET!
	    endlist += "," + str(mytup[0] - 1 + mytup[1])
	    lenlist += "," + str(mytup[1])
	    end_coord = mytup[0] - 1  + mytup[1]

        # remove leading comma
        startlist = re.sub(r',', "", startlist, 1)
        lenlist = re.sub(r',', "", lenlist, 1)
        endlist = re.sub(r',', "", endlist, 1)

	# print line for UCSC bed file
        print (dtran[myisoform][0]  + "\t" + str(start_coord) + "\t" + str(end_coord) + "\t" + newgenename + "\t" + 
	    "1000" + "\t" + dtran[myisoform][2] + "\t" + str(start_coord) + "\t" + str(end_coord) + "\t" + "0" + "\t" + 
	    str(num_blocks) + "\t" + lenlist + "\t" + startlist )

        if verbose:
            print(myisoform)
            print(newgenename)
	    print(startlist)
	    print(lenlist)
            print
