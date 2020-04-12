
# RSDS- A command-line tool for simulating raw RNA-sequencing data 


# Synopsis


# Description
The RNA-seq data simulator described in this study was developed as a command-line interface implemented in Python 3. The tool simulates raw RNA-sequencing data. It takes as input an annotated reference transcript nucleotide sequence file in FASTA format from which reads are simulated. RSDS extracts individual transcript sequences from an index at random in default mode.  The only required arguments are: a user-defined read length, desired number of reads to simulate, library type (single-end or paired-end) and an output file prefix. The simulated RNA-seq reads are written to disk in a FASTQ formatted file with .fastq as file extension. Other parameters to control the properties of the simulated data are available as tuneable settings such as fragment-length distribution, customized Phred-quality score modelling, customized transcript expression profiling, as well as differential transcript expression simulation. 

# Requirements

Python version 3
Numpy
Pyfaidx

# Usage


The options for the program is as follows:

-h --help                           Print the help file and exit
-v --version                        Print the version of the program and exit
-o                                  Prefix for the output files created
-f                                  Reference genome file in FASTA format
-R                                  The desired read length 
