#!/usr/bin/python3
# encoding =UTF-8

import numpy as np
import argparse
from sklearn.mixture import GaussianMixture as GMM
import pickle as pickle
import gzip

parser = argparse.ArgumentParser(description='fragment length distribution modelling')
parser.add_argument('-f', type=str, required=True, help='Input SAM file')
parser.add_argument('-n', type=int, required=True, default=10, help='Number of components')
parser.add_argument('-o', type=str, required=True, help='Output file prefix')

args = parser.parse_args()

samFile = args.f
components = args.n
outfile = args.o


def process_SAM(infile):
    
    FL = []
    f = open(infile, 'r')
    for line in f:
        
        if line[0] != '@':
            parts = line.split('\t')
            TLEN = parts[8]
            FL.append(TLEN)
    
    FL = list(map(int, FL))
    FLS = FL.sort()
    NEG = []
    POS = []
    
    for i in FL:
        if i < 0:
            NEG.append(i)
        else:
            POS.append(i)
    
    data = np.array(POS)
    datalog = np.log(data)
    data2log = datalog.reshape(-1, 1)
    
    return data2log


def percentage(percent, whole):
  return (percent * whole) / 100.0


def optimal_n_components(aic, n):
    n_components = len(aic)


    ll = [i - (index + 1) * 2 for index, i in enumerate(aic)]
    x = [i * percentage(10, n) + j for i, j in zip(ll, range(n_components + 1))]
    minimum = min(x)

    N = [k - minimum for k in x]
    val, idx = min((val, idx) for (idx, val) in enumerate(N))

    return val, idx + 1


def model_fitting(data, n):

    
    ll = [i - (index + 1) * 2 for index, i in enumerate(aic)]
    x = [i * percentage(1, n) + j for i, j in zip(ll, range(n_components + 1))]
    minimum = min(x)
    
    final_data = [k - minimum for k in x]
    val, idx = min((val, idx) for (idx, val) in enumerate(final_data))
    
    return idx + 1


def model_fitting(data, n):
    
    total_obs = len(data)
    aic = []
    n_components_range = range(1, n + 1)
    for n_components in n_components_range:
        gmm = GMM(n_components=n_components, covariance_type='full')
        gmm.fit(data)
        aic.append(gmm.aic(data))

    N = optimal_n_components(aic, len(data))

    
    N = optimal_n_components(aic, total_obs)
    gmm = GMM(n_components=N, covariance_type='full')
    clf = gmm.fit(data)
    mus = clf.means_
    sigmas = clf.covariances_
    weights = clf.weights_
    model = [mus, sigmas, weights]
    
    return model


def main():

    FL_model = model_fitting(process_SAM(samFile), components)
    modelName = outfile + '_p.gzip'
    g = gzip.open(modelName, 'wb')
    pickle.dump(FL_model, g)
    print('Finished!')


if __name__ == '__main__':
    main()




