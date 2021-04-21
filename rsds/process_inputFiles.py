# encoding=utf-8

import pickle as pickle
import numpy as np
import gzip
import itertools
import random
from rsds import distributions, output
from rsds import sequence_handling
import sys
import logging.handlers
from shutil import copyfile


errlog = logging.getLogger("ErrLog")


def defaultfragsize(fragment_size, fragment_std, counts):

	RFS = []
	FS = np.random.normal(fragment_size, fragment_std, 100000).astype(int)
	for i in counts:
		randomFS = random.choices(FS, k=i)
		RFS.append(randomFS)

	return RFS


def find_nearest(array, lower_bound, upper_bound):
	array = np.asarray(array)
	idx_lower = (np.abs(array - lower_bound)).argmin()
	idx_upper = (np.abs(array - upper_bound)).argmin()
	
	return idx_lower, idx_upper


def parseIndexRef(indexFile):
	"""
	Description:
	Read in sequence data from reference index FASTA file returns a list of transcript IDs
	offset, seqLen, position
	Parameters
	 - indexFile (str): The index file generated by the program, written to the current directory
	Return: The function returns a list of tuples containing the transcript id, start and end offset of the transcript sequence
	"""
	ref_inds = []
	filt_ref_inds = []

	try:

		fai = open(indexFile, 'r')
	except BaseException:

		errlog.error('Cannot find indexed reference file. Please provide a reference FASTA file')
		sys.exit('Cannot find indexed reference file. Please provide a reference FASTA file')

	for line in fai:
		splt = line[:-1].split('\t')
		header = '>' + splt[0]
		seqLen = int(splt[1])
		offset = int(splt[2])
		lineLn = int(splt[3])
		nLines = seqLen / lineLn

		if seqLen % lineLn != 0:
			nLines += 1
		ref_inds.append([header, offset, offset + seqLen + nLines, seqLen])
	for i in ref_inds:
		if i[3] >= 400:
			filt_ref_inds.append(i)
	for x in filt_ref_inds:
		x.pop(3)
	fai.close()

	return filt_ref_inds


def default_simulation(refFile, refindex, readtot):

	counts = []
	transcript_sequences = []
	NB_counts = distributions.negative_binomial()
	counts_NB = np.random.choice(NB_counts, size=readtot, replace=True).tolist()
	scaled_counts = sequence_handling.scalereadnum(counts_NB, readtot)
	samptransids = random.choices(refindex, k=len(scaled_counts))

	for i in samptransids:
		p = sequence_handling.processTransIDs(refFile, [i])
		for id, seq in p.items():
			transcript_sequences.append(seq)
	counts.append(scaled_counts)

	return transcript_sequences, counts


def samplingtranscripts(ids, readtot):

	""""
	Description: This function randomly sample from all reference transcripts
	Parameters: ids (list of tuples) It takes as input all reference transcripts offsets
	Returns: This function returns a subset of transcript ids to be sampled from
	"""
	random.seed(seed)
	numreads = readtot
	sampledtranscripts = random.sample(ids, numreads)

	return sampledtranscripts


def proc_FLmodel(file, n):

	with gzip.open(file, 'rb') as f:

		model = pickle.load(f)

		mus = model[0].flatten()
		sigma = np.sqrt(model[1].flatten())
		weights = model[2]
		aic = model[3]
		bic = model[4]
		size = [int(round(i * n)) for i in weights]
		np.random.seed(1234)
		sample = []
		for m, sc, si in zip(mus, sigma, size):
			N = np.random.normal(loc=m, scale=sc, size=si)

			sample.append(N.tolist())

		merged_dist = list(itertools.chain.from_iterable(sample))
		res = np.sort(np.round(np.exp(merged_dist)).astype(int)).tolist()

		return res


def sample_target(dist, readlen, seqlen, counts):

	bounds = find_nearest(dist, readlen, seqlen)
	target_sample = dist[bounds[0]:bounds[1]]
	randomFS = random.choices(target_sample, k=counts)
	
	return randomFS


def process_models(profile):

	transcript_ID = []
	transcript_count = []
	transcript_propcount = []

	for i in profile:

		rec_id = i[0]
		trans_id = rec_id.split('|')
		id = [trans_id[0], i[1], i[2]]
		transcript_ID.append(id)
		rec_counts = i[3]
		transcript_count.append(rec_counts)
		rec_propcount = i[4]
		transcript_propcount.append(rec_propcount)

	return transcript_ID, transcript_count, transcript_propcount


def get_trans_sequences(transcript_offsets, reference):

	transcript_sequences = []
	for i in transcript_offsets:
		p = sequence_handling.processTransIDs(reference, [i])
		for id, seq in p.items():
			transcript_sequences.append(seq)

	return transcript_sequences


def diffmode(model, reference, readtot):

	with open(model, 'rb') as file:
		profile = pickle.load(file)

		background = process_models(profile[0].tolist())
		control = process_models(profile[1].tolist())
		experiment = process_models(profile[2].tolist())

		b_prop = np.rint(np.multiply(background[2], readtot)).astype(int)
		background_counts = [1 if i==0 else i for i in b_prop]
		background_seq = get_trans_sequences(background[0], reference)

		c_prop = np.rint(np.multiply(control[2], readtot)).astype(int)
		control_counts = [1 if i==0 else i for i in c_prop]
		control_seq = get_trans_sequences(control[0], reference)

		e_prop = np.rint(np.multiply(experiment[2], readtot)).astype(int)
		experiment_counts = [1 if i==0 else i for i in e_prop]
		experiment_seq = get_trans_sequences(experiment[0], reference)
		data = {'background':(background_seq, background_counts), 'control':(control_seq, control_counts), 'experiment': (experiment_seq, experiment_counts)}

	return data


def makefraglendist(fmodel, readtot, readlen, counts, transcript_sequences):

	fraglendist = []

	FRS = proc_FLmodel(fmodel, readtot)
	for i, j in zip(transcript_sequences, counts):
		randomFS = sample_target(FRS, readlen, len(i), j)
		fraglendist.append(randomFS)

	return fraglendist


def compilefastqrecord(readlen, reference, refindex, qmodel, filename, kwargs, model=None, diff=None, readtot=None, fragmodel=None):

	global counts

	if model==None and diff==None:

		default_data = default_simulation(reference, refindex, readtot)
		defcounts = list(itertools.chain.from_iterable(default_data[1]))

		if kwargs is 'se':
			sim_data = sequence_handling.assemble_reads(default_data[0], defcounts, readlen, qmodel, kwargs)
			output.write_fastq(filename, sim_data, single_end=True)

		elif kwargs is 'pe':
			if fragmodel is None:
				counts = defaultfragsize(250, 25, defcounts)
			else:
				counts = makefraglendist(fragmodel, readtot, readlen, defcounts, default_data[0])
			sim_data = sequence_handling.assemble_reads(default_data[0], counts, readlen, qmodel, kwargs)
			output.write_fastq('rsdsv0.1', sim_data, single_end=False)

	elif model != None and diff==None:

		with gzip.open(model, 'rb') as file:
			profile = pickle.load(file)
			profile_data = process_models(profile)
			if readtot==None:
				counts = profile[1]
			else:
				counts = np.rint(np.multiply(profile_data[2], readtot)).astype(int)
				counts = [1 if i==0 else i for i in counts]
			seq = get_trans_sequences(profile_data[0], reference)

		if kwargs is 'se':
			sequence_handling.assemble_reads(seq, counts, readlen, qmodel, kwargs)
		elif kwargs is 'pe':
			if fragmodel is None:
				counts = defaultfragsize(250, 25, counts)
			else:
				counts = makefraglendist(fragmodel, readtot, readlen, counts, seq)
			sequence_handling.assemble_reads(seq, counts, readlen, qmodel, kwargs)

	elif model == None and diff != None:

		diff_data = diffmode(diff, reference, readtot)
		allfilenames = []
		if kwargs is 'se':
			seq, counts = diff_data['background']
			background_data = sequence_handling.assemble_reads(seq, counts, readlen, qmodel, kwargs)
			outfilename = output.write_fastq('grp1', background_data, single_end=True)
			cfile = output.generatefilename(name='grp2', zipped=True, single_end=True)
			copyfile(outfilename, cfile)
			allfilenames.append(outfilename)
			allfilenames.append(cfile)
			seq, counts = diff_data['control']
			control_data = sequence_handling.assemble_reads(seq, counts, readlen, qmodel, kwargs)
			seq, counts = diff_data['experiment']
			exp_data = sequence_handling.assemble_reads(seq, counts, readlen, qmodel, kwargs)
			output.add_simreads(allfilenames, control_data, single_end=True, simdata2=exp_data)

		elif kwargs is 'pe':

			seq, counts = diff_data['background']
			if fragmodel is None:
				fcounts = defaultfragsize(250, 25, counts)
			else:
				fcounts = makefraglendist(fragmodel, readtot, readlen, counts, seq)
			background_data = sequence_handling.assemble_reads(seq, fcounts, readlen, qmodel, kwargs)
			outfile1 = output.write_fastq('grp1', background_data, single_end=False)
			cfile = output.generatefilename(name='grp2', zipped=True, single_end=False)
			copyfile(outfile1[0], cfile[0])
			copyfile(outfile1[1], cfile[1])
			allfilenames.append(outfile1)
			allfilenames.append(cfile)
			seq, counts = diff_data['control']
			control_data = sequence_handling.assemble_reads(seq, fcounts, readlen, qmodel, kwargs)
			seq, counts = diff_data['experiment']
			exp_data = sequence_handling.assemble_reads(seq, fcounts, readlen, qmodel, kwargs)
			output.add_simreads(allfilenames, control_data, single_end=False, simdata2=exp_data)




