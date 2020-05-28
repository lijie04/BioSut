#!/usr/bin/env python3

# -*- coding:utf-8 -*-

from distutils.core import setup

import os
from bioutil.system import files

def version():
	vers_file = os.path.dirname(__file__) + '/VERSION'
	files.check_exist(vers_file)

	return open(vers_file).readline().strip()

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
