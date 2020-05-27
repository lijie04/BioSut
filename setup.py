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
	author = 'Jie Li',
	author_email = 'jeveylijie@gmail.com',
	url = 'None',
	description = 'bioutil lib for bio related bioinformatics operations.',
	license = 'MIT3.0',
	platform = 'Linux',
	packages = ['bioutil']
	)
