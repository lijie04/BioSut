"""
The :mod:`biosut.gt_file` includes functions relate to file operations.
"""

# Author: Jie Li (mm.jlli6t@gmail.com)
# License: GNU v3.0

import os
import sys
import gzip
import re
from .gt_path import abs_path


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
		#logger.error('File * %s * does not exists.', f)
		sys.exit('File * %s * does not exists.' % f)

def check_file_exist(*files, check_empty:bool=False):
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
			#logger.error('File * %s * does not exists.', f)
			sys.exit('File * %s * does not exists.' % f)
		if check_empty:
			check_file_empty(f)

def check_file_empty(*files):
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
			#logger.error('File * %s * is empty.', f)
			sys.exit('File * %s * is empty.' % f)

def get_file_prefix(file_in:str=None, times:int=1, split_symbol:str='.', include_path:bool=False):
	"""
	Get prefix of file, e.g. file is test.fa, then return test

	Parameters
	----------
	file_in : str
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
	#check_file_exist(file_in) #check outside, to escape file not exist but still want the prefix
	file_in = abs_path(file_in).split(split_symbol)
	file_in = split_symbol.join(file_in[0:len(file_in)-times])
	if include_path:return file_in
	return os.path.basename(file_in)

def get_seqfile_prefix(seqin:str=None):
	"""
	Get prefix of a sequence file.
	e.g. file is someting_1.fastq.gz, return something.

	Parameter
	---------
	seqin : str
		Input FASTQ/FASTA format file.

	Return
	------
	str:
		Return a string indicate as the prefix of seqfile.
	"""

	seqin = os.path.basename(seqin)
	return re.sub('.\d.fastq.gz|.\d.fq.gz|.\d.fastq|.\d.fq|.\d.fa.gz|.\d.fasta.gz|.\d.fa|.\d.fasta', '', seqin)

# Deprecated, to add times parameter
#	for i in range(times):
#		if include_path:return os.path.splitext(f)[0]
#		return os.path.splitext(os.path.basename(f))[0]

def get_file_path(*files):
	"""
	Get absolute path of input and return.

	Parameters
	----------
	*files : str
		Input file (s)

	Return
	------
	str
		Absolute path (s) of your input file (s).
	"""
	final_paths = [os.path.dirname(abs_path(f)) for f in files]

	if len(final_paths) == 1:return final_paths[0]
	return final_paths

def perfect_open(file_in:str):
	"""
	Make a perfect open for file

	Parameters
	----------
	file_in : str
		input file to open

	Returns
	-------
	str
		Return file handle (s)
	"""
	#final_files = []
	#for f in file_in:
		#final_files.append(gzip.open(f, 'rt')) if '.gz' in f else final_files.append(open(f, 'r'))
	if '.gz' in file_in:return gzip.open(file_in, 'rt')
	return open(file_in, 'r')

def close_file(*file_handle):
	"""
	Close file handles.

	Parameters
	----------
	file_handle : str
		Input file handles (s)

	Results
	-------
		Close all file handles.
	"""
	for f in file_handle:
		f.close()

def find_files(Dir, suffix='fa'):
	"""
	Find files with suffix.

	Parameters
	----------
	Dir : str
		Input directory to find files
	suffix : str, default 'fa'
		suffix of files

	Returns
	-------
	list
		A list of found files with full path.
	"""
	return [Dir+'/'+f for f in os.listdir(Dir) if f.endswith(suffix)]
