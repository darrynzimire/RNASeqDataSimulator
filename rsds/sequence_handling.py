# encoding=utf-8

import numpy as np
import random
from Bio.Seq import Seq
from rsds import output
import SequenceContainer
from rsds.main import se_class, sqmodel


# SE_class = SequenceContainer.ReadContainer(m.readlen, m.sqmodel, m.SE_RATE)


def sample_qualscore(sqmodel):

	(myQual, myErrors) = se_class.getSequencingErrors(sqmodel)

	return myQual


def scalereadnum(read_counts, n):
	sc = []
	scale = []
	total = sum(read_counts)
	for i in read_counts:
		x = i / total
		scale.append(x)
	for i in scale:
		y = n * i
		sc.append(round(y))
	scaled_counts = [1 if x == 0 else x for x in sc]

	return scaled_counts


def getseq(reference, key, start=1, end=None):
	"""
	Description
	Get a sequence by key coordinates are 1-based and end is inclusive
	Parameters:
		key:
		start:
		end:
	Returns:
	"""

	if end != None and end < start:
		return ""
	start -= 1
	seek = start

	# if seek is past sequence then return empty sequence
	if seek >= end:
		return ""

	# seek to beginning
	with open(reference, 'r') as infile:

		infile.seek(seek)

		# read until end of sequence
		header = ''
		seq = []
		if end == None:
			lenNeeded = util.INF
		else:
			lenNeeded = end - start

		len2 = 0
		while len2 < lenNeeded:
			line = infile.readline()
			if line.startswith(">") or len(line) == 0:
				break
			seq.append(header + line.rstrip())
			len2 += len(seq[-1])
			if len2 > lenNeeded:
				seq[-1] = seq[-1][:-int(len2 - lenNeeded)]
				break
		seq = "".join(seq)

		return seq


def processTransIDs(refFile, ids):

	""""
	Description:
	This function take as input a list of transcript ids and converts it to a dictionary
	Parameters:
		ids (list of tuples): List of transcript ids
	Returns: The function returns a dictionary of transcript id as key and start and end position as value
	"""

	Transseq = []
	header = []

	transcriptID = {i: [j, k] for i, j, k in ids}
	ID = transcriptID.keys()
	for k in ID:
		header.append(k)
	pos = transcriptID.values()
	for i in pos:
		start = i[0]
		end = i[1]
		seq = getseq(refFile, ID, start, end)
		Transseq.append(seq)

	new_dict = {k: v for k, v in zip(header, Transseq)}

	return new_dict


def GenerateRead(seq, readLen, n, *args):

	"""
	Description:
	This function truncates transcript sequences by a specified read length.
	Parameters:
	:param seq: Transcript sequence randomly sampled from the input reference transcriptome file
	:param readLen: The user-specified read length

	:return: The function returns a list of all truncated sequences
	"""

	seqLen = len(seq)

	spos = []
	epos = []
	for ag in args:

		if ag == 'SE':

			nmax = seqLen - readLen - 1
			v = np.round(np.random.uniform(low=0, high=nmax, size=n))
			startpos = list(random.choices(v, k=n))
			endpos = [i + readLen for i in startpos]
			spos.append(startpos)
			epos.append(endpos)

		elif ag == 'PE':

			nmax = [seqLen - i - 1 for i in n]
			v = np.round(np.random.uniform(low=0, high=nmax, size=None))
			startpos = list(random.choices(v, k=len(n)))
			endpos = [i + j for i, j in zip(startpos, n)]
			spos.append(startpos)
			epos.append(endpos)

	return spos, epos


def reverse_complement(inputread):
	s = Seq(str(inputread))
	read = s.reverse_complement()
	return read


def process_reads_PE(readlen, fragment):

	R1 = []
	R2 = []
	prob = str(np.random.rand(1)).lstrip('[').rstrip(']')
	read1 = ''.join(map(str, fragment[:readlen]))
	read2 = str(reverse_complement(fragment[-readlen:]))
	if float(prob) < 0.5:
		R1.append(read2)
		R2.append(read1)
	else:
		R1.append(read1)
		R2.append(read2)

	return R1, R2


def get_reads(record):
	start = record[0]
	end = record[1]
	sequence = record[2]
	reads = []
	for s, e in zip(start, end):
		r = sequence[int(s):int(e)]
		reads.append(r)

	return reads


def assemble_reads(Seq, counts, readlen, qmodel, kwargs):

	sim_data = []
	if kwargs == 'se':
		for seq, r in zip(Seq, counts):
			readinfo = GenerateRead(seq, readlen, r, 'SE')
			startpos = readinfo[0]
			endpos = readinfo[1]
			for index, (i, j) in enumerate(zip(startpos[0], endpos[0])):
				header = output.assemble_Illumina_line(instrument='rsdsv0.1', single_end=True)
				read = seq[int(i):int(j)]
				quality_string = sample_qualscore(qmodel)
				records = header, read, quality_string
				sim_data.append(records)

	elif kwargs == 'pe':
		for seq, r in zip(Seq, counts):

			readinfo = GenerateRead(seq, readlen, r, 'PE')
			startpos = readinfo[0]
			endpos = readinfo[1]
			R1 = []
			R2 = []
			for index, (i, j) in enumerate(zip(startpos[0], endpos[0])):
				read = seq[int(i):int(j)]
				data = process_reads_PE(readlen, read)
				R1.append(''.join(data[0]))
				R2.append(''.join(data[1]))
			for index, (r1, r2) in enumerate(zip(R1, R2)):
				header = output.assemble_Illumina_line(instrument='rsdsv0.1', single_end=False)
				quality_string = sample_qualscore(qmodel)
				reads = (r1, r2)
				records = (header, reads, quality_string)
				sim_data.append(records)

	return sim_data

