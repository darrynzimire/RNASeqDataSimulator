# encoding=utf-8


import pickle as pickle
import numpy as np
import gzip
import itertools


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
	res = np.round(np.exp(merged_dist))

	return res


def proc_tx_expmodel(model):

	transcript_ID = []
	transcript_count = []
	transcript_propcount = []
	file = open(model, 'rb')
	profile = list(pickle.load(file, encoding='latin1'))
	profile = [list(i) for i in profile]

	for i in profile:
		transcript_id = tuple(i[0:3])
		transcript_count = i[3]
		transcript_propcount = i[4]
		transcript_ID.append(transcript_id)
		transcript_count.append(transcript_count)
		transcript_propcount.append(transcript_propcount)

	return transcript_ID, transcript_count, transcript_propcount


def proc_qualmodel(model):
	pass


def process_fastq(file):
	pass



def process_countTable(table):
	pass


def process_SAM(file):
	pass


def proc_FASTA():
	pass