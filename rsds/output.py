# encoding=utf-8

"""This module handles all textual output"""

import gzip


def sequence_identifier(readlen, index):

	header = '{}{} {} {}{}'.format('@RSDS_v0.1.', index, index, 'length=', str(readlen))
	return header


def write_fastq(outfilename, header, read, kwargs):

		if kwargs == 'se':
			with open(outfilename + '.fastq', 'w') as handle:
				handle.write('{}\n{}\n+\n'.format(header, read))

		elif kwargs == 'pe':
			with open(outfilename + '_R1.fastq', 'w') as f1, open(outfilename + '_R2.fastq','w') as f2:
				f1.write('{}\n{}\n+'.format(header, read))
				f2.write('{}\n{}\n+'.format(header, read))


def genelist(model):

	genelist = []
	for i in model:
		# print(i)
		line = i[0].split('|')
		# print(line)
		gene_name = line[5]
		genelist.append(gene_name)

	return genelist


def generatefilename():
	pass


def assemblestatsdata():
	pass


def write_ground_truth():
	pass

def simulation_statistics():
	pass


def write_analysis_reportPDF():
	pass


def plotting():
	pass