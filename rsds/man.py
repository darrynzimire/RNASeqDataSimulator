# encoding=utf-8


def manpage():

	text =  'RSDS: version 1.0.0 \n'\
			'author: Darryn Zimire \n'\
			'author email: darrenzim@gmail.com/ darrynzim@sun.ac.za \n'\
			'source: http://github.com/darrynzimire/RNASeqDataSimulator \n'\
			'License: MIT \n'\
		   '\n'\
		   'RNA-Seq Data Simulator (RSDS) is a command-line interface implemented in Python 3. \n' \
		   'The tool simulates raw RNA-sequencing data by emulating characteristics of real RNA-seq data. \n' \
		   'Parameters to control the properties of the simulated data are available as tuneable settings\n' \
		   'such as fragment-length distribution, customized Phred-quality score modelling, customized \n' \
		   'transcript expression profiling, as well as the simulation of paired-end and single-end data.\n'

	print(text)
	'\n'
	print('program: rsds-simulate\n'
	"\n"
	  'This is the main program to simulate RNA-seq data. \n'
	  "\n"
	  'usage: rsds-simulate <[options]> \n'
	  "\n"
	  'Options: \n'
		  '\n'
		  '-r           read length \n'
		  '-f           reference transcriptome FASTA file \n'
		  '-n           sequencing depth (total number of reads to simulate \n'
		  '-se          single-end simulation \n'
		  '-pe          paired-end simulation \n'
		  '-q           quality-score model \n'
		  '-fl          fragment-length distribution model to sample fragment sizes \n'
		  '-o           output files prefix \n'
	  "\n"
	  'program: rsds-learn-qmodel\n'
	  "\n"
	  'This tool creates a position-wise distribution of quality values from a fastq file \n'
	 'It creates a .qmodel.p file which can be passed to RSDS to simulate Phred quality scores. \n'

		"\n"
	  'usage: rsds-learn-qmodel <[options]>\n'
		  "\n"
	  'options: \n'
		  '\n'
		  '-i           input FASTQ file \n'
		  '-o           output file prefix \n'
		  "\n"
	  'program: rsds-learn-FLmodel:\n'
			  "\n"

	'This tool models the fragment-length distribution pattern derived from the aligned insert length.\n'
	'A Gaussian mixture model is fitted to the extracted data and the best-fit model selected is returned.\n'
	'The resulting model can be used in simulation to sample fragment-sizes for paired-end simulation.\n'
			  "\n"
	  'usage: rsds-learn-FLmodel <[options]>\n'
				  "\n"
	  'options:\n'
		  '\n'
		  '-f           input SAM file \n'
		  '-n           number of components (Gaussians to fit) default=10 \n'
		  '-o           output file prefix \n'
				  "\n"
	  'program: rsds-learn-profile:\n'
				  "\n"
	'RSDS provides a feature to reproduce a empirical transcript expression profile.The program takes\n'
	'in a count matrix with transcript identities and expression estimate in the form of counts\n'
	'and returns a model which can be used in simulation to reproduce that particular expression profile \n'
	
	  "\n"
	  'usage: rsds-learn-profile <[options]>\n'
				  "\n"
	  'options:\n'
		  '\n'
		  '-f           reference FASTA file\n'
		  '-c           count matrix in text format\n'
		  '-o           output prefix\n'
			'\n'
	  'program: rsds-diff:\n'
				  "\n"
	  'description:\n'
				  "\n"
	  'usage:\n'
				  "\n"
	  'options:\n'

		)


def main():

	manpage()


if __name__ == '__main__':
	main()