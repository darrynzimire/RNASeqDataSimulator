# encoding=utf-8

import numpy as np
import pickle as pickle
import itertools


def proc_FLmodel(file, n):

	np.random.seed(1234)
	modelled_dist = []
	model = pickle.load(file)
	weights = model[0].flatten()
	means = model[1].flatten()
	covars = model[2].flatten()
	size = round(sum([i * n for i in weights]))

	for m, s in zip(means, covars):
		N = np.round(np.exp(np.random.normal(loc=m, scale=s, size=int(size)))).astype(int)
		modelled_dist.append(N.flatten().tolist())
	merged_dist = list(itertools.chain.from_iterable(modelled_dist))

	return merged_dist


def proc_tx_expmodel(model):

	profile = pickle.load(model)
	ids = profile[0]
	counts = profile[1]
	propcount = profile[2]

	return ids, counts, propcount


def proc_qualmodel():
	pass