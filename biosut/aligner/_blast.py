"""
The :mod:`biosut.aligner._blast` run blast tool box.
"""

__author__ = 'Jie Li'
__copyright__ = 'Copyright 2023'
__credits__ = 'Jie Li'
__license__ = 'GPLv3.0'
__maintainer__ = 'Jie Li'
__email__ = 'jlli6t near gmail.com'

import os

class blast(object):
	"""Wrapper for running blast"""

	def __init__(self, fasta, db, ):
		
