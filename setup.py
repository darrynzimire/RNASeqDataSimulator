# encoding=utf-8

import setuptools

with open('README.md', 'r', encoding='utf-8') as fh:
	long_description = fh.read()

setuptools.setup(
	name="rsds",
	version="1.0.0",
	author="Darryn Zimire",
	author_email="darrynzim@sun.ac.za",
	description="Simulator for RNA-sequencing datasets to inform experimental design.",
	long_description=long_description,
	long_description_content_type='markdown',
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
		"Topic :: Scientific/Engineering :: Bioinformatics",
		],
	install_requires=["numpy", "matplotlib", "seaborn",
					  "pyfaidx", "pandas", "biopython", "scipy",
					  "scikit-learn", 'setuptools'],

	packages=setuptools.find_packages(),
	entry_points={
		'console_scripts': [
			'rsds=rsds.man:main',
			'rsds-simulate = rsds.__main__:main',
			'rsds-learn-qmodel = rsds.tools.learn_Qmodel:main',
			'rsds-learn-FLmodel = rsds.tools.learn_FLmodel:main',
			'rsds-learn-profile = rsds.tools.Transcript_Expression_Profiling:main'

		]},
)

