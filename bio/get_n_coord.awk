#!/bin/awk -f

# About:
# find coordinates of the Ns in a fasta file

# Useage:
# $ cat example/test.fa | get_n_coord.awk

BEGIN{flag=0}{if ($0 ~ />/) {if (flag==1) {print nstart"\t"nend;}; print; pos=0; flag=0} else {for (i=1; i<=length; i++) { pos++; c=substr($0,i,1); if (c=="N" && flag==0) {nstart=pos; nend=pos; flag=1;} else if (c=="N") {nend=pos} else if (c!="N" && flag==1) {print nstart"\t"nend; flag=0}}}}END{if (flag==1) {print nstart"\t"nend;}}
