#!/usr/bin/python3
# encoding =UTF-8

import pickle as pickle

def process_countmodel(model):

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



def process_fraglenmodel(model):
	pass


def process_fastq(file):
	pass


def process_countTable(table):
	pass


def process_SAM(file):
	pass
