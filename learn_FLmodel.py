#!/usr/bin/python3
# encoding =UTF-8

import os
import numpy as np
import argparse
from sklearn.mixture import GaussianMixture as GMM

import pickle as pickle
import gzip

parser = argparse.ArgumentParser(description='fragment length distribution modelling')
parser.add_argument('-f', type=str, required=True, help='Input SAM file')
parser.add_argument('-n', type=int, required=True, help='Number of components')
parser.add_argument('-o', type=str, required=True, help='Output file prefix')

args = parser.parse_args()

samFile = args.f
components = args.n
outfile = args.o


def process_SAM(infile):
    
    samFile = ''
    for file in os.listdir('.'):
        if file.endswith('.sam'):
            samFile = os.path.join('.', file)
    
    FL = []
    f = open(samFile, 'r')
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


def model_fitting(data, n):
    
    global best_gmm
    
    lowest_aic = np.infty
    aic = []
    n_components_range = range(1, n + 1)
    for n_components in n_components_range:
        gmm = GMM(n_components=n_components, covariance_type='full')
        gmm.fit(data)
        aic.append(gmm.aic(data))
        if aic[-1] < lowest_aic:
            lowest_aic = aic[-1]
            best_gmm = gmm
    
    clf = best_gmm
    mus = clf.means_
    sigmas = clf.covariances_
    weights = clf.weights_
    model = [weights, mus, sigmas]
    
    return model


def main():
    
    FL_model = model_fitting(process_SAM(samFile), components)
    modelName = outfile + '_p.gzip'
    g = gzip.open(modelName, 'wb')
    pickle.dump(FL_model, g)


if __name__ == '__main__':
    main()




