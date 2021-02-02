# encoding=utf-8


import pickle as pickle
import numpy as np
import gzip
import itertools
import random


def find_nearest(array, lower_bound, upper_bound):
	array = np.asarray(array)
	idx_lower = (np.abs(array - lower_bound)).argmin()
	idx_upper = (np.abs(array - upper_bound)).argmin()
	
	return idx_lower, idx_upper


def getmax_seqlen():
	pass


def proc_FLmodel(file, n):

	f = gzip.open(file, 'rb')
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
	res = np.sort(np.round(np.exp(merged_dist)))
	
	return res


def sample_target(dist, readlen, seqlen, counts):
	
	bounds = find_nearest(dist, readlen, seqlen)
	target_sample = dist[bounds[0]:bounds[1]]
	randomFS = random.choices(target_sample, k=counts)
	
	return randomFS


def proc_tx_expmodel(model):

	transcript_ID = []
	transcript_count = []
	transcript_propcount = []
	file = open(model, 'rb')
	profile = pickle.load(file, encoding='latin1')

	for i in profile:
		rec_id = i[0]
		trans_id = rec_id.split('|')
		id = (trans_id[0], i[1], i[2])
		transcript_ID.append(id)
		rec_counts = np.rint(i[3]).astype(int)
		transcript_count.append(rec_counts)
		rec_propcount = i[4]
		transcript_propcount.append(rec_propcount)

	return transcript_ID, transcript_count, transcript_propcount


def proc_DEmodel(model):
	
	file = open(model, 'rb')
	model = pickle.load(file)
	group_A = model[0]
	group_B = model[1]
	
	return group_A, group_B


def proc_qualmodel(model):
	pass


def process_fastq(file):
	pass


def process_SAM(file):
	pass


def proc_FASTA():
	pass