#!/usr/bin/python3
# encoding = UTF-8

from setuptools import setup

setup(
	name='RSDS',
	description='Simulator for sequencing depth, read length and base-call quality in RNA-sequencing data',
	long_description=open('README.rst').read(),
	version='0.1',
	author='Darryn Zimire',
	author_email='darrynzim@sun.ac.za',
	license='MIT',
	install_requires= ['numpy', 'pyfaidx', 'scipy', 'pandas'],

)