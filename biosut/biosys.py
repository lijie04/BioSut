#####################################################################
#																	#
#biosys.py - some system operations related to biology?				#
#																	#
#####################################################################

__author__ = 'M.M Jie Li'
__copyright__ = 'Copyright 2018'
__credits__ = ['M.M Jie Li']
__license__ = 'GPL3'
__maintainer__ = ['M.M Jie Li']
__email__ = 'mm.jlli6t near gmail.com'

import os
import sys
import logging
import subprocess as sp

logger = logging.getLogger(__name__)

class path:
	@classmethod
	def check_exist(cls, *paths, check_empty : bool=False):
		"""
		Check if path is exists, exits if path does not exist.

		Parameters
		----------
		paths : str (s)
			Path (s) that will be checked.
		check_empty : bool, default True
			Whether to check path emptiness.

		Return
		-------
			Return full path (s).
		"""
		final_paths = []
		for p in paths:
			p = cls.abs_path(p)
			final_paths.append(p)
			if not os.path.exists(p):
				logger.error('Path *%s* does not exists.', p)
				sys.exit()
			if check_empty:
				cls.check_empty(p)
		if len(final_paths) == 1:return final_paths[0]
		return final_paths

	@classmethod
	def sure_exist(cls, *paths):
		"""
		Check if path exists, path will be created if not exists.

		Parameters
		----------
		paths : str
			Path (s) that will be checked
			if path is not existed, path will be created.

		Results
		--------
			Create and return full path (s).
		"""

		final_paths = []
		for p in paths:
			p = path.abs_path(p)
			final_paths.append(p)
			if not os.path.exists(p):
				try:
					os.makedirs(p)
				except OSError as e:
					logger.error('Path *%s* is not creatable.', p, exc_info=True)
					sys.exit()
		if len(final_paths) == 1:return final_paths[0]
		return final_paths
	
	def check_empty(*dirs):
		"""
		Check if directory(s) is empty.
		
		Parameters
		-----------
		dirs : str
			directory (s) to check.

		Results:
		--------
			Exit and report error msg while input directory (s) is empty.
		"""

		for d in dirs:
			if not os.listdir(d):
				logger.error('Directory %s is empty.', d)
				sys.exit()

	def real_path(*paths):
		"""
		Return real path, link file/path will convert to solid original destination path.

		Parameters
		----------
		paths : str
			Input path (s).

		Return
		------
			Real path (s).

		Examples
		--------
		>>> from biosut.biosys import path
		>>> a = './../../bucket'
		>>> b = '/usr/bin/python'
		>>> c = 'aaa -> ../../aaa' # this is a link file.
		>>> final_paths = path.real_path(a, b)
		>>> print(final_paths)
		['/full/path/to/./../../bucket', '/usr/bin/python']
		>>> final_paths = path.real_path(b)
		>>> print(final_paths)
		'/usr/bin/python'
		>>> final_paths = path.real_path(c)
		'/full/path/to/../../aaa'
		"""

		final_paths = [os.path.realpath(p) for p in paths]

		if len(final_paths) == 1:
			return final_paths[0]
		else:
			return final_paths
	
	def abs_path(*paths):
		"""
		Return absulote path (s), link file/path will keep as linked destination path.
		The only difference between real_path and abs_path is about the linked file (s)/path (s).

		Parameters
		----------
		paths : str
			Input path (s).

		Return
		------
		str
			abs path (s)

		Examples
		--------
		>>> from biosut.biosys import path
		>>> a = './../../bucket'
		>>> b = '/usr/bin/python'
		>>> c = 'aaa -> ../../aaa' # this is a link file.
		>>> final_paths = path.abs_path(a, b)
		>>> print(final_paths)
		['/full/path/to/./../../bucket', '/usr/bin/python']
		>>> final_paths = path.abs_path(b)
		>>> print(final_paths)
		'/usr/bin/python'
		>>> final_paths = path.abs_path(c)
		'/full/path/to/aaa'
		"""
		final_paths = [os.path.abspath(p) for p in paths]
		if len(final_paths) == 1:
			return final_paths[0]
		else:
			return final_paths

	@classmethod
	def get_path(cls, *files):
		"""
		Get absolute path of input and return.
		
		Parameters
		----------
		paths : str
			Input path (s)

		Return
		-------
		str
			Absolute path (s) of your input file (s).
		"""
		
		final_paths = [os.path.dirname(cls.abs_path(f)) for f in files]
		
		if len(final_paths) == 1:
			return final_paths[0]
		else:
			return final_paths

	def check_program(*prog):
		"""
		Check whether program (s) exists in system path.

		Paramters
		---------
		prog : str
			Program (s) that will be check.
		
		Result
		------
			Exit if program (s) is not exist in system path.
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

		Parameters
		----------
		db_v :	str
			Database name that will be check

		Result
		------
			Exit process if db is not exist.
		"""
		db = os.environ.get(db_v)
		if db:
			cls.check_empty(db)
		else:
			logger.error('Did not find %s'%db)
			sys.exit()
		return db

	def exe_proc(cmd, shell : bool=True):
		"""
		Executing your command.
		
		Parameters
		----------
		cmd:str, or list
			Command will be executed.
		shell : bool, default True
			Set to False if cmd input is a list.
		
		Results
		-------
			Execute command and return output result and error messages.
		"""
		proc = sp.Popen(cmd, shell=shell, stdout=sp.PIPE, stderr=sp.PIPE)
		out, err = proc.communicate()
		if proc.returncode != 0:
			logger.error('Error encountered while executing:\n%s\nError message:\n%s\n' %(cmd, err))
			sys.exit()
		return out, err

### class files
class files:
	
	## deprecated, redundanted.
	def _check(f):
		"""
		Check if file exists, if file is not exists, program exit.

		Parameter
		---------
		f : str
			Input file that will be check.
		"""

		if not os.path.isfile(f):
			logger.error('File * %s * does not exists.', f)
			sys.exit()

	@classmethod
	def check_exist(cls, *fls, check_empty : bool=False):
		"""
		Check if file (s) exists, exit if file not existed.

		Parameters
		----------
		fls : str
			Input file (s) to check existance.
		check_empty : bool, default False.
			File emptiness will be checked.

		Result
		------
			Process will be killed if file is not exist.
			If `check_empty` is True, then process will be killed if file is empty.
		"""
		for	f in fls:
			if not os.path.isfile(f):
				logger.error('File * %s * does not exists.', f)
				sys.exit()
			if check_empty:
				cls.check_empty(f)

	def check_empty(*fls):
		"""
		Check if file is empty.

		Parameters
		----------
		fls : str
			Input file (s) to check emptiness.

		Results
		-------
			Exit and report errors if file is empty.
		"""
		for f in fls:	
			if not os.path.getsize(f):
				logger.error('File * %s * is empty.', f)
				sys.exit()

	def get_prefix(f, times : int=1, include_path : bool=False):
		"""
		Get prefix of file, e.g. file is test.fa, then return test
		
		Parameters
		----------
		f : str
			Input file, relative path or abusolute path.
		times : int, default 1.
			How many times supposed to chop str behind symbol '.'
		include_path : bool, default False
			Whether to include path or not.

		Returns
		-------
			Return prefix.
		"""
		f = path.abs_path(f)
		for i in range(times):
			if include_path:
				return os.path.splitext(f)[0]
			else:
				return os.path.splitext(os.path.basename(f))[0]

	def perfect_open(f):
		"""
		Make a perfect open for file

		Returns
		-------
		str
			Return file handle
		"""
		
		if '.gz' in f:
			import gzip
			return gzip.open(f, 'rt')
		else:
			return open(f, 'r')



