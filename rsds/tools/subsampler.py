# encoding=utf-8

import numpy as np
import random
import argparse
import gzip
import sys

# sampling from paired-end files
# writing output to a folder
# checks to make sure sampling is truely random
# using filename function to write file names
# writing out a log file
# Perhaps an option for sampling reads through a theoretical distribution


parser = argparse.ArgumentParser(description='FASTQ subsampler')
parser.add_argument('--f1',         type=str, required=False, metavar='<str>', help='* input_read1.fq (.gz)"')
parser.add_argument('--f2',         type=str, required=False, metavar='<str>', help='* input_read2.fq (.gz)"')
parser.add_argument('--depth',      type=int, required=False, metavar='<str>', default=5000, help='')
parser.add_argument('--samples',    type=int, required=False, metavar='<str>', default=10, help='')
parser.add_argument('--libtype',    type=str, required=False, metavar='<str>', default='paired', help='')
parser.add_argument('--seed',       type=int, required=False, metavar='<str>', default=12345, help='')
parser.add_argument('--output',     type=str, required=False, metavar='<str>', default='rsdsv1.0.samples', help='')


args = parser.parse_args()
f1 = args.f1
f2 = args.f2
depth = args.depth
samples = args.samples
libtype = args.libtype
seed = args.seed
output = args.output


def randind(min, max, n):

	dist = np.round(np.random.uniform(low=min, high=max, size=max))
	sample = np.random.choice(dist, size=n, replace=False)
	rand_ind = np.array(sorted(sample)).astype(int)
	print(type(rand_ind))
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


print("counting records....")
with gzip.open(f1) as input:
	num_lines = sum([1 for line in input])
total_records = int(num_lines / 4)
#
# if args.fraction:
# 	args.number = int(total_records * args.fraction)

print("sampling " + str(depth) + " out of " + str(total_records) + " records")

output_files = []
output_sequence_sets = []
for i in range(samples):
	output_files.append(open(output + "." + str(i) + '.fastq', "w"))
	output_sequence_sets.append(set(random.sample(range(total_records + 1), depth)))

record_number = 0
with gzip.open(f1) as input:
		for line1 in input:
			line2 = next(input).decode()
			line3 = next(input).decode()
			line4 = next(input).decode()
			for i, output in enumerate(output_files):
				if record_number in output_sequence_sets[i]:
						output.write(line1.decode())
						output.write(line2)
						output.write(line3)
						output.write(line4)
				record_number += 1
				if record_number % 100000 == 0:
					print(str((int(record_number / total_records)) * 100) + " % done")


for output in output_files:
	output.close()
print("done!")





