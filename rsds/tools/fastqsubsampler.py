# encoding=utf-8

import gzip
import random
import argparse
import os
import sys
import time

directory = './RSDSv1.0_samples'
file_path = os.path.join(directory)
if not os.path.isdir(directory):
	os.mkdir(directory)

parser = argparse.ArgumentParser()
parser.add_argument('-i1', type=str, required=True,     help="name of FASTQ file 1")
parser.add_argument('-i2', type=str, required=False,    help="name of FASTQ file 2")
parser.add_argument('-s',  type=int, required=False,    help="number of samples", default=3)
parser.add_argument('-n',  type=int, required=False,    help="total amount of reads to sample", default=100)
parser.add_argument('-l',  type=str, required=False,    help="library type: paired or single", default='pe')


args = parser.parse_args()

f1 = args.i1
f2 = args.i2
samples = args.s
number = args.n
libtype = args.l


def generatefilename(name, zipped, n, single_end):

	filenames  = []
	file_path = os.path.join(directory, name)
	for i in range(n):

		basename = "{}_{}".format(file_path, 's' + str(i+1))
		lib_type1 = basename + "_R1.fastq" + (".gz" if zipped else "")
		if single_end:
			filenames.append(lib_type1)
		else:
			lib_type2 = basename + "_R2.fastq" + (".gz" if zipped else "")
			filenames.append((lib_type1, lib_type2))

	return filenames


def getreadseq(outfile, rand_records, total_records,  file1, file2=None, single_end=False):

	if single_end is True:

		with gzip.open(file1) as input:
				suba = gzip.open(outfile, 'wt')
				rec_no = - 1
				for r in rand_records:
					while rec_no < r:
						rec_no += 1
						for i in range(4): input.readline().decode()
					for i in range(4):
						suba.write(input.readline().decode())
					rec_no += 1
					if rec_no % 10 == 0:
						print(str((int(rec_no / total_records)) * 100) + " % done")

	elif single_end is False and file2 is not None:
		with gzip.open(file1) as fha, gzip.open(file2) as fhb:
			suba, subb = gzip.open(outfile[0], 'wt'), gzip.open(outfile[1], 'wt')
			rec_no = - 1
			for rr in rand_records:
				while rec_no < rr:
					rec_no += 1
					for i in range(4): fha.readline().decode()
					for i in range(4): fhb.readline().decode()
				for i in range(4):
					suba.write(fha.readline().decode())
					subb.write(fhb.readline().decode())
				rec_no += 1


def write_records(single_end):

	print("counting records....")

	records = sum(1 for _ in gzip.open(f1, 'rb')) // 4
	print("sampling " + str(number) + " out of " + str(records) + " records")

	if single_end:

		output_files = generatefilename(name='rsds_v1.0', zipped=True, n=samples, single_end=True)
		for i, f in zip((range(samples)), output_files):
			rand_records = tuple(sorted(set(random.sample(range(records + 1), number))))
			rr = {f: rand_records}
			for key, value in rr.items():
				getreadseq(key, value, records, f1, f2, single_end=True)
	elif not single_end:
		counter = 0
		output_files = generatefilename(name='rsds_v1.0', zipped=True, n=samples, single_end=False)
		for i, f in zip((range(samples)), output_files):
			rand_records = tuple(sorted(set(random.sample(range(records + 1), number))))
			rr = {f: rand_records}
			for key, value in rr.items():
				getreadseq(key, value, records, f1, f2, single_end=False)
				counter += 1
		print(counter)


def main():

	if libtype == 'se':
		write_records(single_end=True)
	elif libtype == 'pe':
		write_records(single_end=False)


if __name__ == '__main__':
	main()
