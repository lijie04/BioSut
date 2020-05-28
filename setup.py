#!/usr/bin/env python3

# -*- coding:utf-8 -*-

from distutils.core import setup

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


if __name__ == '__main__':

	setup(
		name = 'bioutil',
		version = version(),
		author = 'Jie Li',
		author_email = 'jeveylijie@gmail.com',
		url = 'Not now',
		description = 'bioutil lib for bio related bioinformatics operations.',
		license = 'Not now',
		platform = 'Linux',
		packages = ['bioutil']
		)
