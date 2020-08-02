# -*- coding:utf-8 -*-

from setuptools import setup, find_packages
import os

def version():
	vers1 = os.path.dirname(os.path.realpath(__file__)) + '/VERSION'
	vers2 = os.path.dirname(os.path.realpath(__file__)) + './biosut/VERSION'
	if os.path.isfile(vers1):return open(vers1).readline().strip()
	if os.path.isfile(vers2):return open(vers2).readline().strip()
	return 'unknown version'


setup(
	name = 'biosut',
	description = 'biology suite for bioinformatics operations.',
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
