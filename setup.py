# -*- coding:utf-8 -*-

from setuptools import setup, find_packages

import os
import sys

def version():
	vers_file = os.path.dirname(os.path.realpath(__file__)) + '/VERSION'
	
	if os.path.isfile(vers_file):
		version = open(vers_file).readline().strip()

		print('You are installing bioutil %s' % version)
		return version
	else:
		sys.exit('Dont know your package version, didnt find VERSION file')


setup(
	name = 'biosut',
	version = version(),
	description = 'bioutil library for bio related bioinformatics operations.',
	url='https://github.com/jlli6t/biosut',
	author = 'M.M Jie Li',
	author_email = 'mm.jlli6t@gmail.com',

	classifiers = [
					'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
					'Programming Language :: Python :: 3 :: Only',
					'Operating System :: Unix',
					],
	packages = find_packages(),
	python_requires = '>=3.6',

	)


