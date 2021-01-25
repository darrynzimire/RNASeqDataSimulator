# encoding=utf-8

import os
import gzip
import numpy as np
import argparse
import pickle as pickle
from rsds.probability import DiscreteDistribution


def parseFQ(inf):
    print('reading ' + inf + '...')
    if inf[-3:] == '.gz':
        print('detected gzip suffix...')
        f = gzip.open(inf, 'r')
    else:
        f = open(inf, 'r')
    
    IS_SAM = False
    if inf[-4:] == '.sam':
        print('detected sam input...')
        IS_SAM = True
    
    rRead = 0
    actual_readlen = 0
    qDict = {}
    while True:
        
        if IS_SAM:
            data4 = f.readline()
            if not len(data4):
                break
            try:
                data4 = data4.split('\t')[10]
            except IndexError:
                break
        # need to add some input checking here? Yup, probably.
        else:
            data1 = f.readline()
            data2 = f.readline()
            data3 = f.readline()
            data4 = f.readline()
            if not all([data1, data2, data3, data4]):
                break
        
        if actual_readlen == 0:
            if inf[-3:] != '.gz' and not IS_SAM:
                totalSize = os.path.getsize(inf)
                entrySize = sum([len(n) for n in [data1, data2, data3, data4]])
                print('estimated number of reads in file:', int(float(totalSize) / entrySize))
            actual_readlen = len(data4) - 1
            print('assuming read length is uniform...')
            print('detected read length (from first read found):', actual_readlen)
            priorQ = np.zeros([actual_readlen, RQ])
            totalQ = [None] + [np.zeros([RQ, RQ]) for n in range(actual_readlen - 1)]
        
        # sanity-check readlengths
        if len(data4) - 1 != actual_readlen:
            print('skipping read with unexpected length...')
            continue
        
        for i in range(len(data4) - 1):
            q = ord(data4[i]) - offQ
            qDict[q] = True
            if i == 0:
                priorQ[i][q] += 1
            else:
                totalQ[i][prevQ, q] += 1
                priorQ[i][q] += 1
            prevQ = q
        
        rRead += 1
        if rRead % PRINT_EVERY == 0:
            print("\r Analysing {:>10} reads".format(rRead), end="")
        
        if MAX_READS > 0 and rRead >= MAX_READS:
            break
    f.close()
    
    print('computing probabilities...')
    probQ = [None] + [[[0. for m in range(RQ)] for n in range(RQ)] for p in range(actual_readlen - 1)]
    for p in range(1, actual_readlen):
        for i in range(RQ):
            rowSum = float(np.sum(totalQ[p][i, :])) + PROB_SMOOTH * RQ
            if rowSum <= 0.:
                continue
            for j in range(RQ):
                probQ[p][i][j] = (totalQ[p][i][j] + PROB_SMOOTH) / rowSum
    
    initQ = [[0. for m in range(RQ)] for n in range(actual_readlen)]
    for i in range(actual_readlen):
        rowSum = float(np.sum(priorQ[i, :])) + INIT_SMOOTH * RQ
        if rowSum <= 0.:
            continue
        for j in range(RQ):
            initQ[i][j] = (priorQ[i][j] + INIT_SMOOTH) / rowSum
    
    print('estimating average error rate via simulation...')
    Qscores = list(range(RQ))
    
    initDistByPos = [DiscreteDistribution(initQ[i], Qscores) for i in range(len(initQ))]
    probDistByPosByPrevQ = [None]
    for i in range(1, len(initQ)):
        probDistByPosByPrevQ.append([])
        for j in range(len(initQ[0])):
            if np.sum(probQ[i][j]) <= 0.:  # if we don't have sufficient data for a transition, use the previous qscore
                probDistByPosByPrevQ[-1].append(DiscreteDistribution([1], [Qscores[j]], degenerateVal=Qscores[j]))
            else:
                probDistByPosByPrevQ[-1].append(DiscreteDistribution(probQ[i][j], Qscores))
    
    countDict = {}
    for q in Qscores:
        countDict[q] = 0
    for samp in range(1, N_SAMP + 1):
        if samp % PRINT_EVERY == 0:
           print("\r Analysing {:>10} reads ".format(samp), end="")
        "\n"
        myQ = initDistByPos[0].sample()
        countDict[myQ] += 1
        for i in range(1, len(initQ)):
            myQ = probDistByPosByPrevQ[i][myQ].sample()
            countDict[myQ] += 1
    
    totBases = float(sum(countDict.values()))
    avgError = 0.
    for k in sorted(countDict.keys()):
        eVal = 10. ** (-k / 10.)
        
        avgError += eVal * (countDict[k] / totBases)
    "\n"
    print('\nAVG ERROR RATE:', avgError)
    
    return (initQ, probQ, avgError)


parser = argparse.ArgumentParser(description='learn-qmodel.py')
parser.add_argument('-i', type=str, required=True, metavar='<str>', help="* input_read1.fq (.gz)")
parser.add_argument('-o', type=str, required=True, metavar='<str>', help="* output.p")
parser.add_argument('-i2', type=str, required=False, metavar='<str>', default=None,
                    help="input_read2.fq (.gz) / input_read2.sam")
parser.add_argument('-p', type=str, required=False, metavar='<str>', default=None, help="input_alignment.pileup")
parser.add_argument('-q', type=int, required=False, metavar='<int>', default=33, help="quality score offset [33]")
parser.add_argument('-Q', type=int, required=False, metavar='<int>', default=41, help="maximum quality score [41]")
parser.add_argument('-n', type=int, required=False, metavar='<int>', default=-1,
                    help="maximum number of reads to process [all]")
parser.add_argument('-s', type=int, required=False, metavar='<int>', default=1000000,
                    help="number of simulation iterations [1000000]")
args = parser.parse_args()

(INF, OUF, offQ, maxQ, MAX_READS, N_SAMP) = (args.i, args.o, args.q, args.Q, args.n, args.s)
(INF2, PILEUP) = (args.i2, args.p)

RQ = maxQ + 1

INIT_SMOOTH = 0.
PROB_SMOOTH = 0.
PRINT_EVERY = 10000


def main():
    Qscores = list(range(RQ))
    if INF2 == None:
        (initQ, probQ, avgError) = parseFQ(INF)
    else:
        (initQ, probQ, avgError1) = parseFQ(INF)
        (initQ2, probQ2, avgError2) = parseFQ(INF2)
        avgError = (avgError1 + avgError2) / 2.
    
    if PILEUP == None:
        
        print('Using default sequencing error parameters...')
        
        # sequencing substitution transition probabilities
        SSE_PROB = [[0., 0.4918, 0.3377, 0.1705],
                    [0.5238, 0., 0.2661, 0.2101],
                    [0.3754, 0.2355, 0., 0.3890],
                    [0.2505, 0.2552, 0.4942, 0.]]
        # if a sequencing error occurs, what are the odds it's an indel?
        SIE_RATE = 0.01
        # sequencing indel error length distribution
        SIE_PROB = [0.999, 0.001]
        SIE_VAL = [1, 2]
        # if a sequencing indel error occurs, what are the odds it's an insertion as opposed to a deletion?
        SIE_INS_FREQ = 0.4
        # if a sequencing insertion error occurs, what's the probability of it being an A, C, G, T...
        SIE_INS_NUCL = [0.25, 0.25, 0.25, 0.25]
    
    else:
        print('\nPileup parsing coming soon!\n')
        exit(1)
    
    errorParams = [SSE_PROB, SIE_RATE, SIE_PROB, SIE_VAL, SIE_INS_FREQ, SIE_INS_NUCL]
    
    print('saving model...')
    if INF2 == None:
        pickle.dump([initQ, probQ, Qscores, offQ, avgError, errorParams], open(OUF + '.p', 'wb'))
    else:
        pickle.dump([initQ, probQ, initQ2, probQ2, Qscores, offQ, avgError, errorParams], open(OUF + '.p', 'wb'))


if __name__ == '__main__':
    main()
