
#####################################################################
#																	#
#system.py - system operations related functions				#
#																	#
#####################################################################

__author__ = 'Irene Jie Li'
__copyright__ = 'Copyright 2018'
__credits__ = ['Irene Jie Li']
#__license__ = 'GPL3'  choose the proper open license
__maintainer__ = ['Irene Jie Li']
__email__ = 'irene.jie.li6 at gmail.com'

import os
import sys
import subprocess
import logging

logger = logging.getLogger(__name__)

class path:

	@classmethod
	def check_exist(cls, *paths, check_empty=False):
		"""
		Check if path is exists,
		program exits if path does not exist.

		Parameters:
		p: Path taht will be checked
		check: bool, check path empty or not.

		Result:
		This will return all full paths.
		"""
		new_paths = []
		for p in paths:
			p = os.path.realpath(p)
			new_paths.append(p)
			if not os.path.exists(p):
				logger.error('Path *%s* does not exists.', p)
				sys.exit()
		
			if check_empty:
				cls.check_empty(p)
		if len(new_paths) == 1:return new_paths[0]
		return new_paths


	def sure_exist(*paths):
		"""
		Check if path is exists, path will be created if it
		if not exists then create path.

		Parameters:
		p: Path that will be sure exists

		Result:
		This will return all full paths.
		"""

		new_paths = []
		for p in paths:
			p = os.path.realpath(p)
			new_paths.append(p)
			if not os.path.exists(p):
				try:
					os.makedirs(p)
				except OSError as e:
					logger.error('Path *%s* is not creatable.', p, exc_info=True)
					sys.exit()
		if len(new_paths) == 1:return new_paths[0]
		return new_paths
		
	def check_empty(*dirs):
		"""
		Check if directory is empty.
		
		Parameters:
		d: directory to check

		results:
		If it is empty, report errors.
		"""
		for d in dirs:
			if not os.listdir(d):
				logger.error('Directory %s is empty.', d)
				sys.exit()


	def realpath(p):
		return os.path.realpath(p)


	def check_program(*prog):
		"""
		Check if program exists or not.

		Paramters:
		prog: program that will be check
		"""
		for p in prog:
			code = subprocess.run(['which', p], stdout=subprocess.PIPE, stderr=subprocess.STDOUT).returncode

			if code:
				logger.error('Program * %s * is not found', p)
				sys.exit()

	@classmethod
	def check_db(cls, db_v):
		"""
		Check if db exists or not.

		Parameters:
		db:	db that will be check

		Result:
		if db not exist, program quit.
		"""
		db = os.environ.get(db_v)
		if db:
			cls.check_empty(db)
		else:
			logger.error('Did not find %s'%db)
			sys.exit()
		return db

#	deprecated, use os.path.realpath instead
#	def make_full_path(self):
#		#print(path)
#		if self.path.startswith('/'):
#			self.path = self.path
#		elif self.path.starswith('.'):
#			self.path = os.path.join(os.getcwd(), self.path)
#		else:
#			self.path = os.path.join(os.getcwd(), self.path)
#		self.make_sure_path_exists(self.path)
#		return self.path

class files:
	
	## deprecated, redundanted.
	def _check(f):
		"""
		Check if file exists, if file is not exists, program exit.

		Parameter:
		f: file that will be check.
		"""

		if not os.path.isfile(f):
			logger.error('File * %s * does not exists.', f)
			sys.exit()

#	@classmethod
#	def check_exist(cls, in_file, check_empty=False):
		"""
		Check if file exist or empty.

		Parameters:
		in_file:	file that will be checked, can be a file or file list
		check_empty:	bool, check if file empty or not.

		Results:
		If file is not exists, or empty, report errors and exit.

		"""

#		if type(in_file) == str:
#			if not os.path.isfile(in_file):
##				logger.error('File *%s* does not exists', in_file)
#				sys.exit()
#			if check_empty:
#				files.check_empty(in_file)
#		else:
#			if type(in_file) == list:
#				for f in in_file:
#					check_file_exist(f)
	
	@classmethod
	def check_exist(cls, *kargs, check_empty=False):
		for	k in kargs:
			if not os.path.isfile(k):
				logger.error('File * %s * does not exists.', k)
				sys.exit()
			if check_empty:
				cls.check_empty(k)
		

	def check_empty(f):
		"""
		Check if file is empty.

		Parameters:
		f: the file that will be checked.

		Results:
		If file is empty, report errors.
		"""
	
		if not os.path.getsize(f):
			logger.error('File * %s * is empty.', f)
			sys.exit()



