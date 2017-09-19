###bio
(some rubbish in here :)

- **aaseq_2_aakmer.pl** produce kmers from an amino acid sequence; translate kmers into nucleotides if specified (dependencies: BioPerl)

- **break_200_N.sh** a script to take a fasta file and break it into smaller entries on stretches of "N" characters (default: 200 N's). 

- **exonconcat.cpp** a simple program that, given a sequence and a file of ranges (say, representing exons), extracts the substring of ranges from the parent sequence and concatenates them

- **fastajoinlines** for each entry in a fasta file, join the sequence portion of the file onto a single line (removing all but one newline character)

- **get_n_coord.awk** given a fasta file, get the coordinates of "N" characters or stretches thereof

- **get_nt.cpp** a simple program to get substrings from a sequence file

- **get_perc_overlap.pl** from a standard tab-delimited blast output file, get the *cumulative* subject and query coverages (dependencies: bedtools)

- **get_top_blast_hit.pl** given a blast file with degenerate query IDs, output a blast file of uniq query IDs (by taking the first one seen)	

- **gtf2ucsc_bed.py** a script which converts exon ranges in gtf format into UCSC appropriate bed format (as specified [here](https://genome.ucsc.edu/FAQ/FAQformat.html#format1))

- **humangencode.txt** a text file of codon to amino acid mappings

- **nt2aa.cpp** a simple program to convert nucleotides into amino acids given the codon mappings

- **plot_data.r** given a file of data in columnar form, plot it using R 

- **partition_genome.pl** partition a file of genomic ranges into smaller subdivisions

- **read_sam_flag.pl** a script to display the FLAG field of a SAM (Sequence Alignment/Map format) file in human-readable form

- **revc.cpp** a simple program to reverse complement a nucleotide sequence
