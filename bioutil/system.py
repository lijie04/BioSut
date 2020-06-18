
#####################################################################
#																	#
#system.py - system operations related functions				#
#																	#
#####################################################################

__author__ = 'Jie Li'
__copyright__ = 'Copyright 2018'
__credits__ = ['Jie Li']
#__license__ = 'GPL3'  choose the proper open license
__maintainer__ = ['Jie Li']
__email__ = 'jeveylijie at gmail.com'

import os
import sys
import logging
import subprocess as sp

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
			p = path.abs_path(p)
			new_paths.append(p)
			if not os.path.exists(p):
				logger.error('Path *%s* does not exists.', p)
				sys.exit()
		
			if check_empty:
				cls.check_empty(p)
		if len(new_paths) == 1:return new_paths[0]
		return new_paths

	@classmethod
	def sure_exist(cls, *paths):
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
			p = path.abs_path(p)
			new_paths.append(p)
			if not os.path.exists(p):
				try:
					os.makedirs(p)
				except OSError as e:
					logger.error('Path *%s* is not creatable.', p, exc_info=True)
					sys.exit()
		if len(new_paths) == 1:return new_paths[0]
		return new_paths
	
	def exe_proc(cmd, shell=True):
		"""
		Executing your command.
		
		Parameters:
		-----------
		cmd:str/list
			input command
		shell=bool
			True, set to False if cmd is a list.
		
		Results:
		Return output result and error messages.
		"""

		proc = sp.Popen(cmd, shell=shell, stdout=sp.PIPE, stderr=sp.PIPE)
		out, err = proc.communicate()
		if proc.returncode != 0:
			logger.error('Error encountered while executing:\n%s\nError message:\n%s\n' %(cmd, err))
			sys.exit()
		return out, err


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


	def real_path(*paths):
		"""
		This will return real path.

		Parameters:
		-----------
		paths:str
			input path (s).

		Return:
		-------
		str
			real path of input path
		"""

		new_paths = [os.path.realpath(p) for p in paths]

		if len(new_paths) == 1:
			return new_paths[0]
		else:
			return new_paths
	
	def abs_path(*paths):
		"""
		This will return abusolute path.

		Parameters:
		-----------
		paths:str
			input path (s).

		Return:
		-------
		str
			abs path of input path.
		"""
		
		new_paths = [os.path.abspath(p) for p in paths]
		if len(new_paths) == 1:
			return new_paths[0]
		else:
			return new_paths


	@classmethod
	def get_path(cls, *files):
		"""
		Get absolute path of your file and return it.
		
		Parameters:
		-----------
		paths:str
			input path (s)

		Return:
		-------
		str
			Absolute path of your input file
		"""
		
		final_paths = [os.path.dirname(cls.abs_path(f)) for f in files]
		
		if len(final_paths) == 1:
			return final_paths[0]
		else:
			return final_paths


	def check_program(*prog):
		"""
		Check if program exists or not.

		Paramters:
		prog: program that will be check
		"""
		for p in prog:
			code = sp.run(['which', p], stdout=sp.PIPE, stderr=sp.STDOUT).returncode

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


### class files
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


	def get_prefix(f, times=1, include_path=False):
		"""
		Get prefix of file, e.g. file is test.fa, then return test
		
		Parameters:
		-----------
		f:str
			input file, relative path or abusolute path.
		times:int
			how many times I supposed to chop str behind symbol '.'
		include_path:bool
			bool value to include absolute path or not.

		Returns:
		--------
			return prefix(s).
		"""
		f = path.abs_path(f)
		for i in range(times):
			if include_path:
				f = os.path.splitext(f)[0]
			else:
				f = os.path.splitext(os.path.basename(f))[0]
		return f

	def perfect_open(f):
		"""
		Make a perfect open for file

		Returns
		-------
		str
			Return a file handle
		"""
		
		if '.gz' in f:
			import gzip
			return gzip.open(f, 'rt')
		else:
			return open(f, 'r')

