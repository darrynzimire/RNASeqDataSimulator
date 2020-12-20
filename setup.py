# encoding = UTF-8

import setuptools

with open('README.md', 'r', encoding='utf-8') as fh:
	long_description=fh.read()

setuptools.setup(
	name="rsds",
	version="0.0.1",
	author="Darryn Zimire",
	author_email="darrynzim@sun.ac.za",
	description="Simulator for RNA-seq datasets for experimental design.",
	long_description=long_description,
	long_description_content_type='text/markdown',
	license="MIT",
	url='https://github.com/darrenzimire/RNASeqDataSimulator',
	classifiers=[
		"Intended Audience :: Science/Research",
		"License :: MIT License",
		"Natural Language :: English",
		"Programming Language :: Python :: 3",
		"Programming Language :: Python :: 3.5",
		"Programming Language :: Python :: 3.6",
		"Programming Language :: Python :: 3.7",
		"Programming Language :: Python :: 3.8",
		"Topic :: Scientific/Engineering :: Bio-Informatics",
		],
	install_requires=["numpy", "matplotlib", "seaborn",
					  "pyfaidx", "pandas", "biopython", "scipy",
					  "scikit-learn"],
	packages=setuptools.find_packages()

	)

