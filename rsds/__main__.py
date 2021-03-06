# encoding=utf-8

import os, warnings
import sys
import pyfaidx
from rsds import SequenceContainer
from rsds import process_inputFiles
from rsds import distributions
from rsds.utilities import sequence_handling
import argparse
import logging.handlers
from datetime import datetime
from rsds import man, output


if not sys.warnoptions:
	warnings.simplefilter("default")  # Change the filter in this process
	os.environ["PYTHONWARNINGS"] = "default"  # Also affect subprocesses
	warnings.simplefilter("ignore", ResourceWarning)

# Make a global logging object.
errlog = logging.getLogger("ErrLog")

# Set logging level, and write everything to a file
errlog.setLevel(logging.DEBUG)
LOG_FILENAME = './err.log'
h = logging.FileHandler(LOG_FILENAME, 'w')
f = logging.Formatter("%(levelname)s %(asctime)s %(funcName)s %(lineno)d %(message)s")
h.setFormatter(f)
errlog.addHandler(h)


def get_arguments():
	
	parser = argparse.ArgumentParser()
	
	parser.add_argument('-r', 	type=int, 			required=False, 	default=101)
	parser.add_argument('-n', 	type=int, 			required=False)
	parser.add_argument('-f', 	type=str, 			required=False)
	parser.add_argument('-s', 	type=int, 			required=False, 	default=1223)
	parser.add_argument('-o', 	type=str, 			required=False)
	parser.add_argument('-q', 	type=str, 			required=False)
	parser.add_argument('-c', 	type=str, 			required=False)
	parser.add_argument('-er',	type=float, 		required=False, 	default=-1)
	parser.add_argument('-fl',	nargs=2, 		 						default=(250, 25))
	parser.add_argument('-flm', type=str, 			required=False)
	parser.add_argument('-se', action='store_true', required=False)
	parser.add_argument('-pe', action='store_true', required=False)
	parser.add_argument('-de', 	type=str, 			required=False)
	
	return parser


argparser = get_arguments()
args = argparser.parse_args()
(fragment_size, fragment_std) = args.fl
fragmodel = args.flm
ref = args.f
readlen = args.r
readtot = args.n
seed = args.s
outfilename = args.o
sqmodel = args.q
countModel = args.c
diffmodel = args.de
SE_RATE = args.er

start_time = datetime.now()


def main():

	if ref == None:
		man.manpage()
		sys.exit()
	
	else:

		# infile = open(ref, 'r')
		# infile.seek(seek)
		
		errlog.info(print('reading reference file: ' + str(ref) + "\n"))
		errlog.info(print('Indexing reference file....' + "\n"))
		
		basename = str(os.path.basename(outfilename))
		os.symlink(ref, basename)
		pyfaidx.Faidx(basename)
		cwd = os.getcwd()
		indexFile = ''
		for file in os.listdir(cwd):
			if file.endswith('.fai'):
				indexFile = (os.path.join('.', file))
	ref_index = process_inputFiles.parseIndexRef(indexFile)
	SE_CLASS = SequenceContainer.ReadContainer(readlen, sqmodel, SE_RATE)

	if args.se and countModel == None:
		print('running default mode')
		process_inputFiles.compilefastqrecord(readlen, ref, ref_index, outfilename, 'se', model=None, readtot=readtot, fragmodel=None)

	elif args.se and countModel != None:
		process_inputFiles.compilefastqrecord(readlen, ref, ref_index, outfilename,'se', model=countModel, readtot=readtot, fragmodel=None)

	elif args.pe and countModel == None and fragmodel != None:
		process_inputFiles.compilefastqrecord(readlen, ref, ref_index, outfilename, 'pe', model=None, readtot=readtot, fragmodel=fragmodel)

	elif args.pe and countModel != None and fragmodel ==None:
		process_inputFiles.compilefastqrecord(readlen, ref, ref_index, outfilename, 'pe', model=countModel, readtot=readtot, fragmodel=None)

	elif args.pe and countModel != None and fragmodel != None:
		process_inputFiles.compilefastqrecord(readlen, ref, ref_index, outfilename, 'pe', model=countModel, readtot=readtot, fragmodel=fragmodel)

	os.remove(basename)
	os.remove(indexFile)

	errlog.info(print('Simulation is complete'))
	end_time = datetime.now()
	print('Duration: {}'.format(end_time - start_time))


if __name__ == '__main__':
	main()


# 	sample_trans_ids = []
	# 	COUNTS = []
	# 	ID = []
	# 	Seq = []
	#
	# 	if countModel == None:
	# 		profile = process_inputFiles.default_simulation(ref_transcript_ids, readtot)
	# 		errlog.info(print('Simulating single-end reads....' + "\n"))
	# 		errlog.info(print('No transcript profile model detected!!' + "\n"))
	# 		errlog.info(print('Simulating default transcript profile' + "\n"))
	#
	# 		sample_trans_ids.append(profile[0])
	# 		COUNTS.append(profile[1])
	#
	# 	# NB_counts = distributions.negative_binomial()
	# 		# counts_NB = np.random.choice(NB_counts, size=readtot, replace=True).tolist()
	# 		# scaled_counts = sequence_handling.scalereadnum(counts_NB, readtot)
	# 		# samptransids = random.choices(ref_transcript_ids, k=len(scaled_counts))
	# 		# sample_trans_ids.append(samptransids)
	# 		# COUNTS.append(scaled_counts)
	#
	# 	elif countModel != None and readtot == None:
	# 		errlog.info(print('Simulating single-end reads....' + "\n"))
	# 		errlog.info(print('Detected transcript profile model.....' + "\n"))
	# 		errlog.info(print('Simulating empirical transcript profile' + "\n"))
	#
	# 		profile = process_inputFiles.process_models(countModel)
	# 		sample_trans_ids.append(profile[0])
	# 		COUNTS.append(profile[1])
	#
	# 	elif countModel != None and readtot != None:
	# 		profile = process_inputFiles.process_models(countModel)
	# 		counts_s = np.rint(np.multiply(profile[2], readtot)).astype(int)
	# 		COUNTS.append(counts_s)
	# 		sample_trans_ids.append(profile[0])
	#
	# 	if diffmodel != None:
	#
	# 		model = process_inputFiles.proc_DEmodel(diffmodel)
	# 		background = process_inputFiles.process_models(model[0])
	# 		control = process_inputFiles.process_models(model[1])
	# 		experiment = process_inputFiles.process_models(model[2])
	#
	# 	for j in sample_trans_ids:
	# 		p = sequence_handling.processTransIDs(infile, j)
	# 		for id, seq in p.items():
	# 			ID.append(id)
	# 			Seq.append(seq)
	# 	with gzip.open(outfilename + '.fastq.gz', 'wb') as handle:
	# 		for seq, r in zip(Seq, COUNTS[0]):
	#
	# 			readinfo = sequence_handling.GenerateRead(seq, readlen, r, 'SE')
	# 			startpos = readinfo[0]
	# 			endpos = readinfo[1]
	#
	# 			for index, (i, j) in enumerate(zip(startpos[0], endpos[0])):
	# 				header = output.sequence_identifier(readlen, index)
	# 				read = seq[int(i):int(j)]
	# 				q = sequence_handling.sample_qualscore(SE_CLASS, sequencingModel=sqmodel)
	# 				handle.write('{}\n{}\n+\n{}\n'.format(header, read, q).encode())
	#
	# elif args.pe:
	#
	# 	sample_trans_ids = []
	# 	RFS = []
	# 	COUNTS_P = []
	# 	ID = []
	# 	Seq = []
	# 	R1 = []
	# 	R2 = []
	#
	# 	if countModel == None:
	# 		errlog.info(print('Generating paired-end reads.....' + "\n"))
	# 		errlog.info(print('Sampling counts from negative binomial model' + "\n"))
	#
	# 		NB_counts = distributions.negative_binomial()
	# 		counts_NB = np.random.choice(NB_counts, size=readtot, replace=True).tolist()
	# 		counts_p = sequence_handling.scalereadnum(counts_NB, readtot)
	# 		COUNTS_P.append(counts_p)
	# 		sample_trans_ids.append(random.choices(ref_transcript_ids, k=len(COUNTS_P[0])))
	#
	# 	elif countModel != None and readtot == None:
	# 		errlog.info(print('Generating paired-end reads' + "\n"))
	# 		errlog.info(print('Simulating empirical transcript profile.....' + "\n"))
	# 		profile = process_inputFiles.process_models(countModel)
	# 		COUNTS_P.append(profile[1])
	# 		sample_trans_ids.append(profile[0])
	#
	# 	elif countModel != None and readtot != None:
	# 		errlog.info(print('Generating paired-end reads' + "\n"))
	# 		errlog.info(print('Simulating empirical transcript profile.....' + "\n"))
	#
	# 		profile = process_inputFiles.process_models(countModel)
	# 		counts_p = np.rint(np.multiply(profile[2], readtot)).astype(int)
	# 		COUNTS_P.append(counts_p)
	# 		sample_trans_ids.append(profile[0])
	#
	# 	if diffmodel != None:
	# 		model = process_inputFiles.proc_DEmodel(diffmodel)
	# 		background = model[0]
	# 		control = model[1]
	# 		experiment = model[2]
	#
	# 	data = list(itertools.chain.from_iterable(sample_trans_ids))
	# 	for j in data:
	# 		p = sequence_handling.processTransIDs(infile, [j])
	# 		for id, seq in p.items():
	# 			ID.append(id)
	# 			Seq.append(seq)
	#
	# 	if fragmodel != None:
	# 		Fraglen_dist = process_inputFiles.proc_FLmodel(fragmodel, readtot).astype(int)
	#
	# 		for i, j in zip(Seq, COUNTS_P[0]):
	# 			randomFS = process_inputFiles.sample_target(Fraglen_dist, readlen, len(i), j)
	#
	# 			RFS.append(randomFS)
	# 	else:
	# 		FS = np.random.normal(fragment_size, fragment_std, 100000).astype(int).tolist()
	# 		for i in COUNTS_P[0]:
	# 			randomFS = random.choices(FS, k=i)
	# 			RFS.append(randomFS)
	#
	# 	for seq, r in zip(Seq, RFS):
	# 		readinfo = sequence_handling.GenerateRead(seq, r, len(r), 'PE')
	# 		startpos = readinfo[0]
	# 		endpos = readinfo[1]
	# 		for index, (i, j) in enumerate(zip(startpos[0], endpos[0])):
	# 			read = seq[int(i):int(j)]
	# 			data = sequence_handling.process_reads_PE(readlen, read)
	# 			R1.append(''.join(data[0]))
	# 			R2.append(''.join(data[1]))
	# 	with gzip.open(outfilename + '_R1.fastq.gz', 'wb') as f1, gzip.open(outfilename + '_R2.fastq.gz', 'wb') as f2:
	# 		for index, (i, j) in enumerate(zip(R1, R2)):
	# 			id = output.sequence_identifier(readlen, index)
	# 			q1 = sequence_handling.sample_qualscore(SE_CLASS, sequencingModel=sqmodel)
	# 			f1.write('{}\n{}\n+\n{}\n'.format(id, i, q1).encode())
	# 			q2 = sequence_handling.sample_qualscore(SE_CLASS, sequencingModel=sqmodel)
	# 			f2.write('{}\n{}\n+\n{}\n'.format(id, j, q2).encode())
	#
# 	os.remove(basename)
# 	os.remove(indexFile)
#
# 	errlog.info(print('Simulation is complete'))
# 	end_time = datetime.now()
# 	print('Duration: {}'.format(end_time - start_time))
#
#
# if __name__ == '__main__':
# 	main()

