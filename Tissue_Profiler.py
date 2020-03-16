#!/usr/bin/python3
# encoding = UTF-8

import pickle as pickle
import gzip
import argparse
"""
This program takes in a read count expression table and a file output prefix. The table is parsed and puts out a 
"""
parser = argparse.ArgumentParser(description='Tissue_expression_profiler')
parser.add_argument('-f',           type=str,        required=True,        metavar='<str>',       help='Count table')
parser.add_argument('-o',           type=str,        required=True,        metavar='<str>',       help='output file prefix')

args = parser.parse_args()


file = args.f
modelName = args.o


def processCounts(file):

    f = open(file, 'r')
    transcript_ids = []
    Counts = []
    FPKM = []
    header = f.readline()
    for line in f:
        record = line.split('\t')
        # print(record)
        ID = record[1]
        fpkm = record[-1]
        trans_ids = ID[:15]
        counts = round(float(record[4]))
        FPKM.append(fpkm)

        # transcript_ids.append(trans_ids)
        Counts.append(counts)
    print(FPKM)
    Counts = [1 if x == 0 else x for x in Counts]

    # output = modelName + '_p.gzip'
    output = modelName + '.p'

    g = open(output, 'wb')
    pickle.dump(transcript_ids, g)
    pickle.dump(Counts, g)
    g.close()
    f.close()

def main():
    processCounts(file)


if __name__ == '__main__':
    main()
