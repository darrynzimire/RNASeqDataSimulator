# encoding=utf-8

"""This module handles all textual output"""

import platform, socket, re, uuid, json, psutil, logging
import datetime
import random
import gzip
import os
import sys


directory = './RSDSv1.0.sim'
logfile = "RSDSv1.0_simreport.txt"
file_path = os.path.join(directory, logfile)
if not os.path.isdir(directory):
	os.mkdir(directory)


def assemble_Illumina_line(instrument, single_end):

	Illumina_line = "@{}:{}:{}:{}:{}:{}:{} {}:{}:{}" + " "
	run = random.randint(0, 10000)
	flowcell_id = random.randint(0, 10000)
	lane = random.randint(0, 10000)
	tile = random.randint(0, 10000)
	xpos = random.randint(0, 1000000)
	ypos = random.randint(0, 1000000)
	header_1 = Illumina_line.format(instrument, run, flowcell_id, lane, tile, xpos, ypos, 1, 'N', 0)
	if single_end:
		return header_1
	else:
		header_2 = Illumina_line.format(instrument, run, flowcell_id, lane, tile, xpos, ypos, 2, 'N', 0)
		return header_1, header_2


def sequence_identifier(readlen, index):

	header = '{}{} {} {}{}'.format('@RSDS_v1.0.', index, index, 'length=', str(readlen))
	return header


def generatefilename(name, zipped, single_end):

	now = datetime.datetime.now()
	date_time_obj = datetime.datetime.strftime(now, '%Y%m%d-%H%M')
	file_path = os.path.join(directory, name)
	if name is None:
		name = 'RSDSdataset'

	basename = "{}_{}".format(file_path, date_time_obj)
	lib_type1 = basename + "_1.fastq" + (".gz" if zipped else "")
	if single_end:
		return lib_type1
	else:
		lib_type2 = basename + "_2.fastq" + (".gz" if zipped else "")
		return lib_type1, lib_type2


def write_fastq(filename, read_data, single_end):

	global outfilename

	if single_end:
		outfilename = generatefilename(name=filename, zipped=True, single_end=True)
		# file_path = os.path.join(directory, outfilename)
		with open(outfilename, 'a+') as handle:
			for i in read_data:
				handle.write('{}\n{}\n+\n{}\n'.format(i[0], i[1], i[2]))

	if not single_end:
		outfilename = generatefilename(name=filename, zipped=True, single_end=False)
		with open(outfilename[0], 'a+') as f1, open(outfilename[1], 'a') as f2:
			for i in read_data:
				f1.write('{}\n{}\n+\n{}\n'.format(i[0][0], i[1][0], i[2]))
				f2.write('{}\n{}\n+\n{}\n'.format(i[0][1], i[1][1], i[2]))

	return outfilename


def add_simreads(filenames, simdata, single_end, simdata2=None):

	if single_end:
		with open(filenames[0], 'a') as f1, open(filenames[1], 'a') as f2:
			for rec in simdata:
				header = rec[0]
				reads = rec[1]
				qual = rec[2]
				f1.write('{}\n{}\n+\n{}\n'.format(header, reads, qual))
			for rec in simdata2:
				header = rec[0]
				reads = rec[1]
				qual = rec[2]
				f2.write('{}\n{}\n+\n{}\n'.format(header, reads, qual))
	else:
		with open(filenames[0][0], 'a+') as f1, open(filenames[0][1], 'a+') as f2, open(filenames[1][0], 'a+') as f3, open(filenames[1][1], 'a+') as f4:
			for rec in simdata:
				header = rec[0]
				reads = rec[1]
				qual = rec[2]
				f1.write('{}\n{}\n+\n{}\n'.format(header[0], reads[0], qual))
				f2.write('{}\n{}\n+\n{}\n'.format(header[1], reads[1], qual))
			for rec in simdata2:
				header = rec[0]
				reads = rec[1]
				qual = rec[2]
				f3.write('{}\n{}\n+\n{}\n'.format(header[0], reads[0], qual))
				f4.write('{}\n{}\n+\n{}\n'.format(header[1], reads[1], qual))


class Stats(dict):

	def __getattr__(self, attribute):
		return self.get(attribute)
	__setattr__ = dict.__setitem__
	__delattr__ = dict.__delitem__


def collect_simstats(stats):
	pass


def genelist(model):

	genelist = []

	for i in model:
		line = i[0].split('|')
		gene_name = line[5]
		genelist.append(gene_name)

	return genelist


def write_simreport(filepath):

	""" Create and write a report for a run of RSDS"""

	with open(filepath, 'w') as handle:

		now = datetime.datetime.now()
		date_time_obj = datetime.datetime.strftime(now, '%Y-%m-%d %H:%M:%S.%f')
		headerinfo = 'RSDS version 1.0.0 \n' \
					 'Date: ' + date_time_obj[:-7] + '\n'\
					 'Licence: MIT' + '\n'\
					 'Author: Darryn Zimire' + '\n'\
					 'Author email: darrynzim@sun.ac.za\n'\
					 'source: http://github.com/darrynzimire/RNASeqDataSimulator \n'
		description = '\nRSDS is a simulator for RNA-sequencing datasets. It generates synthetic'\
		                ' RNA-sequencing data guided by user-defined parameters.\n'
		'\n'
		parameter_values = '\nSimulation parameter values: \n' \
            '\nrun mode: \n' \
            'read length: \n' \
			'sequencing depth: \n'\
			'library type: \n'\
			'quality profile: \n'
		simulation_results = '\nSimulation results: \n'\


		handle.write(headerinfo)
		handle.write(description)
		handle.write(parameter_values)


write_simreport(file_path)

def getSystemInfo():
    try:
        info={}
        info['platform']=platform.system()
        info['platform-release']=platform.release()
        info['platform-version']=platform.version()
        info['architecture']=platform.machine()
        info['hostname']=socket.gethostname()
        info['ip-address']=socket.gethostbyname(socket.gethostname())
        info['mac-address']=':'.join(re.findall('..', '%012x' % uuid.getnode()))
        info['processor']=platform.processor()
        info['ram']=str(round(psutil.virtual_memory().total / (1024.0 ** 3)))+" GB"
        return json.dumps(info)
    except Exception as e:
        logging.exception(e)


json.loads(getSystemInfo())


def simulation_statistics():
	pass


def write_analysis_reportPDF():
	pass


