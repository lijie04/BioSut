#!/usr/bin/env python3

# -*- coding:utf-8 -*-

from distutils.core import setup

<<<<<<< HEAD
import os

def version():
	bdir = os.path.dirname(os.path.realpath(__file__))
	version = bdir + '/bioutil/VERSION'
	return open(version).readline().strip()

setup(
	name = 'bioutil',
	version = version(),
=======
setup(
	name = 'bioutil',
<<<<<<< HEAD
	version = '0.0.2',
=======
	version = '0.0.1',
>>>>>>> 8627ced38923a418ba1059c7c045700a5f13d292
>>>>>>> 26706aa660cfa7a0dd6f866ff838922b46634e92
	author = 'Jie Li',
	author_email = 'jeveylijie@gmail.com',
	url = 'None',
	description = 'bioutil lib for bio related bioinformatics operations.',
	license = 'MIT3.0',
	platform = 'Linux',
	packages = ['bioutil']
	)
