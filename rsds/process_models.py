# encoding=utf-8

import numpy as np
import itertools
import pickle as pickle
import itertools
import gzip


def proc_FLmodel(file, n):

	f = gzip.open(file, 'r')
	np.random.seed(1234)
	model = list(itertools.chain.from_iterable(pickle.load(f)))
	print(model)
	mus = model[0]
	# print(type(mus))
	# sigma = np.sqrt(model[1]).flatten()
	# weights = model[2]
	# size = [int(round(i * n)) for i in weights]
	# sample = []
	#
	# for m, sc, si in zip(mus, sigma, size):
	# 	N = np.random.normal(loc=m, scale=sc, size=si)
	# 	sample.append(N.tolist())
	#
	# merged_dist = list(map(int, itertools.chain.from_iterable(sample)))
	# merged_dist = np.array(merged_dist)
	# res = np.round(np.exp(merged_dist))
	#
	# return res

def proc_tx_expmodel(model):

	profile = pickle.load(model)
	ids = profile[0]
	counts = profile[1]
	propcount = profile[2]

	return ids, counts, propcount


def proc_qualmodel():
	pass