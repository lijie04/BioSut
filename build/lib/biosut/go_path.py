"""
The :mod:`biosut.go_path` includes functions relate to path operations.
"""

# Author Jie Li (mm.jlli6t@gmail.com)
# License: GNU v3.0

import os
import sys

def check_exist(*paths, check_dir_empty:bool=False):
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
		p = abs_path(p)
		final_paths.append(p)
		if not os.path.exists(p):
			logger.error('Path *%s* does not exists.', p)
			sys.exit()
		if check_empty:
			check_empty(p)
	if len(final_paths) == 1:return final_paths[0]
	return final_paths

def sure_exist(*paths):
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

	Results
	-------
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

	if len(final_paths) == 1:return final_paths[0]
	return final_paths

def abs_path(*paths):
	"""
	Return absulote path (s), link file/path will keep as linked destination path.
	The only different between real_path and abs_path is about the way they dealing with soft links.

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
	if len(final_paths) == 1:return final_paths[0]
	return final_paths

def find_db_path(db_v:str):
	"""
	Check if database is exists or not.

	Parameters
	----------
	db_v :	str
		Database name that will be check

	Returns
	-------
		Return database path If found database and not empty,
		otherwise, program will exits.
	"""

	db = os.environ.get(db_v)
	if db:
		cls.check_empty(db)
        return db
	else:
	    logger.error('Did not find %s'%db)
		sys.exit()
