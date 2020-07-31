"""
The :mod:`biosut.go_file` includes functions relate to file operations.
"""

# Author Jie Li (mm.jlli6t@gmail.com)
# License: GNU v3.0

import os
import sys
import gzip

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

def check_exist(*files, check_file_empty:bool=False):
	"""
	Check if file (s) exists, exit if file not existed.

	Parameters
	----------
	fls : str
		Input file (s) to check existance.
	check_file_empty : bool, default False.
		File emptiness will be checked.

	Result
	------
		Process will be killed if file is not exist.
		If `check_empty` is True, then process will be killed if file is empty.
	"""
	for	f in files:
		if not os.path.isfile(f):
			logger.error('File * %s * does not exists.', f)
			sys.exit()
		if check_empty:
			check_file_empty(f)

def check_empty(*files):
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
	for f in files:
		if not os.path.getsize(f):
			logger.error('File * %s * is empty.', f)
			sys.exit()

def get_prefix(f:str, times:int=1, split_symbol='.', include_path:bool=False):
	"""
	Get prefix of file, e.g. file is test.fa, then return test

	Parameters
	----------
	f : str
		Input file, relative path or abusolute path.
	times : int, default 1.
		How many times supposed to chop str
    split_symbol : str, default '.'
        symbol to use to split file names
	include_path : bool, default False
		Whether to include path or not.

	Returns
	-------
		Return prefix.
	"""
	f = os.path.abs_path(f).split(split_symbol)
    f = split_symbol.join(f[0:len(f)-times-1])
    if include_path:return f
    return os.path.basename(f)

# Deprecated, to add times parameter
#	for i in range(times):
#		if include_path:return os.path.splitext(f)[0]
#		return os.path.splitext(os.path.basename(f))[0]

def get_path(*files):
    """
	Get absolute path of input and return.

	Parameters
	----------
	*files : str
		Input file (s)

	Return
	-------
	str
		Absolute path (s) of your input file (s).
	"""

	final_paths = [os.path.dirname(cls.abs_path(f)) for f in files]

	if len(final_paths) == 1:return final_paths[0]
	return final_paths

def perfect_open(file_in:str):
	"""
	Make a perfect open for file

	Returns
	-------
	str
		Return file handle
	"""

	if '.gz' in f:return gzip.open(f, 'rt')
	return open(f, 'r')
