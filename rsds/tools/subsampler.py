# encoding=utf-8

import numpy as np
import pandas as pd
from Bio import SeqIO
import random
import argparse
import gzip
import io


parser = argparse.ArgumentParser(description='FASTQ subsampler')
parser.add_argument('--f1',         type=str, required=False, metavar='<str>', help='* input_read1.fq (.gz)"')
parser.add_argument('--f2',         type=str, required=False, metavar='<str>', help='* input_read2.fq (.gz)"')
parser.add_argument('--depth',      type=str, required=False, metavar='<str>', default=5000, help='')
parser.add_argument('--libtype',    type=str, required=False, metavar='<str>', default='paired', help='')
parser.add_argument('--seed',       type=str, required=False, metavar='<str>', default=12345, help='')
parser.add_argument('--output',     type=str, required=False, metavar='<str>', default='rsdsv1.0.samples', help='')


args = parser.parse_args()
f1 = args.f1
f2 = args.f2
depth = args.depth
libtype = args.libtype
seed = args.seed
output = args.output


def randind(min, max, n):

	dist = np.round(np.random.uniform(low=min, high=max, size=max))
	sample = np.random.choice(dist, size=n, replace=False)
	rand_ind = sorted(sample)
	print(rand_ind)
	return rand_ind


def generatefilename(name, zipped, n, single_end):

	filenames  = []
	for i in range(n):

		basename = "{}_{}".format(name, 's' + str(i+1))
		lib_type1 = basename + "_1.fastq" + (".gz" if zipped else "")
		if single_end:
			filenames.append(lib_type1)
		else:
			lib_type2 = basename + "_2.fastq" + (".gz" if zipped else "")
			filenames.append((lib_type1, lib_type2))

	return filenames

#
# def getseqrecord(file, index):
#
# 	if file.endswith('.gz'):
# 		f = gzip.open(file, 'rb')
# 	else:
# 		f = open(file, 'r')



print(io.DEFAULT_BUFFER_SIZE)