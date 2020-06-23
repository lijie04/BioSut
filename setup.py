# -*- coding:utf-8 -*-

from setuptools import setup, find_packages
import os

def version():
	vers_file = os.path.dirname(os.path.realpath(__file__)) + '/VERSION'
	
	if os.path.isfile(vers_file):
		return open(vers_file).readline().strip()
	else:
		return none

setup(
	name = 'biosut',
	description = 'biology small suite for related bioinformatics operations.',
	version = version(),
	url='https://github.com/jlli6t/biosut', # optional

	author = 'M.M Jie Li', 
	author_email = 'mm.jlli6t@gmail.com',
	maintainer = 'M.M Jie Li',
	maintainer_email = 'mm.jlli6t@gmail.com',

	classifiers = [
				'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
				'Programming Language :: Python :: 3 :: Only',
				'Operating System :: Unix',
		],
	keywords = 'biology bioinformatics',
	
	packages = find_packages(),
	python_requires = '>=3.6',
)


