# -*- coding:utf-8 -*-

from setuptools import setup, find_packages
from biosut.version import Version

setup(
	name = 'biosut',
	description = 'biology suite for bioinformatics operations.',
	version = Version(),
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
