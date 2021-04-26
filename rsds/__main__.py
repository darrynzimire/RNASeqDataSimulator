# encoding=utf-8

import os, warnings
import sys
import pyfaidx
from rsds import SequenceContainer
from rsds import process_inputFiles
import argparse
import logging.handlers
from datetime import datetime
from rsds import man


if not sys.warnoptions:
	warnings.simplefilter("default")  # Change the filter in this process
	os.environ["PYTHONWARNINGS"] = "default"  # Also affect subprocesses
	warnings.simplefilter("ignore", ResourceWarning)

# Make a global logging object.
errlog = logging.getLogger("ErrLog")

# Set logging level, and write everything to a file

errlog.setLevel(logging.DEBUG)
LOG_FILENAME = './err.log'
h = logging.FileHandler(LOG_FILENAME, 'w')
f = logging.Formatter("%(levelname)s %(asctime)s %(funcName)s %(lineno)d %(message)s")
h.setFormatter(f)
errlog.addHandler(h)


def get_arguments():
	
	parser = argparse.ArgumentParser()
	
	parser.add_argument('-r', 	type=int, 			required=False, 	default=101)
	parser.add_argument('-n', 	type=int, 			required=False)
	parser.add_argument('-f', 	type=str, 			required=False)
	parser.add_argument('-s', 	type=int, 			required=False, 	default=1223)
	parser.add_argument('-o', 	type=str, 			required=False)
	parser.add_argument('-q', 	type=str, 			required=False)
	parser.add_argument('-c', 	type=str, 			required=False)
	parser.add_argument('-er',	type=float, 		required=False, 	default=-1)
	parser.add_argument('-fl',	nargs=2, 		 						default=(250, 25))
	parser.add_argument('-flm', type=str, 			required=False)
	parser.add_argument('-diff', type=str,          required=False)
	parser.add_argument('-se', action='store_true', required=False)
	parser.add_argument('-pe', action='store_true', required=False)

	return parser


argparser = get_arguments()
args = argparser.parse_args()
(fragment_size, fragment_std) = args.fl
fragmodel = args.flm
ref = args.f
readlen = args.r
readtot = args.n
seed = args.s
outfilename = args.o
sqmodel = args.q
countModel = args.c
diffmodel = args.diff
SE_RATE = args.er

start_time = datetime.now()
se_class = SequenceContainer.ReadContainer(readlen, sqmodel, SE_RATE)


def sample_qualscore(sequencingModel):

	(myQual, myErrors) = se_class.getSequencingErrors(sequencingModel)

	return myQual


def main():

	if ref == None:
		man.manpage()
		sys.exit()
	
	else:

		errlog.info(print('reading reference file: ' + str(ref) + "\n"))
		errlog.info(print('Indexing reference file....' + "\n"))
		
		basename = str(os.path.basename(outfilename))
		os.symlink(ref, basename)
		pyfaidx.Faidx(basename)
		cwd = os.getcwd()
		indexFile = ''
		for file in os.listdir(cwd):
			if file.endswith('.fai'):
				indexFile = (os.path.join('.', file))
	ref_index = process_inputFiles.parseIndexRef(indexFile)

	if args.se and countModel == None and diffmodel==None:
		print('running single-end in default mode')
		process_inputFiles.compilefastqrecord(readlen, ref, ref_index, sqmodel, outfilename, 'se', model=None, diff=None, readtot=readtot, fragmodel=None)
	elif args.se and countModel != None and diffmodel==None:
		print('running single-end in empirical profile mode')
		process_inputFiles.compilefastqrecord(readlen, ref, ref_index, sqmodel, outfilename, 'se', model=countModel, diff=None, readtot=readtot, fragmodel=None)
	elif args.se and countModel == None and diffmodel != None:
		print('running single-end in differential expression mode')
		process_inputFiles.compilefastqrecord(readlen, ref, ref_index, sqmodel, outfilename, 'se', model=None,diff=diffmodel, readtot=readtot, fragmodel=None)
	elif args.pe and countModel == None and diffmodel==None and fragmodel == None:
		print('running paired-end in default mode with no FL model')
		process_inputFiles.compilefastqrecord(readlen, ref, ref_index, sqmodel, outfilename, 'pe', model=None, diff=None,  readtot=readtot, fragmodel=None)
	elif args.pe and countModel == None and diffmodel==None and fragmodel != None:
		print('running paired-end in default mode with FL model')
		process_inputFiles.compilefastqrecord(readlen, ref, ref_index, sqmodel, outfilename, 'pe', model=None, diff=None,  readtot=readtot, fragmodel=fragmodel)
	elif args.pe and countModel != None and diffmodel==None and fragmodel ==None:
		print('running paired-end in empirical profile mode with no FL model')
		process_inputFiles.compilefastqrecord(readlen, ref, ref_index, sqmodel, outfilename, 'pe', model=countModel, diff=None,  readtot=readtot, fragmodel=None)
	elif args.pe and countModel != None and diffmodel==None and fragmodel != None:
		print('running paired-end in empirical profile mode with FL model')
		process_inputFiles.compilefastqrecord(readlen, ref, ref_index, sqmodel, outfilename, 'pe', model=None, diff=None,  readtot=readtot, fragmodel=fragmodel)
	elif args.pe and countModel == None and diffmodel !=None and fragmodel ==None:
		print('running paired-end in differential mode with no FL model')
		process_inputFiles.compilefastqrecord(readlen, ref, ref_index, sqmodel, outfilename, 'pe', model=None, diff=diffmodel, readtot=readtot, fragmodel=None)
	elif args.pe and countModel == None and diffmodel !=None and fragmodel != None:
		print('running paired-end in differential mode with FL model')
		process_inputFiles.compilefastqrecord(readlen, ref, ref_index, sqmodel, outfilename, 'pe', model=countModel, diff=diffmodel,  readtot=readtot, fragmodel=fragmodel)

	os.remove(basename)
	os.remove(indexFile)
	errlog.info(print('Simulation is complete'))
	end_time = datetime.now()
	print('Duration: {}'.format(end_time - start_time))


if __name__ == '__main__':
	main()


