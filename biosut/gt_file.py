"""
The :mod:`biosut.gt_file` includes functions relate to file operations.
"""

# Author: Jie Li (mm.jlli6t@gmail.com)
# License: GNU v3.0
# Copyrigth: 2015 -

import os
import sys
import gzip
import re
from .gt_path import abs_path
from loguru import logger

## deprecated, redundanted.
def _check(f):
	"""
	Check if file exists, if file is not exists, program exit.

	Parameter
	---------
	f : str
		Input file that will be check.
	"""

	if not os.path.isfile(f):return False
	return True

def check_file_exist(*files_in, check_empty:bool=False):
	"""
	Check if file (s) exists, exit if file not existed.

	Parameters
	----------
	files_in : str
		Input file (s) to check existance.
	check_empty : bool, default False.
		File emptiness will be checked.

	Result
	------
		return a boolean value.
	"""
	for	fl in files_in:
		if not os.path.isfile(fl):
			logger.error(f"{fl} is not exist.")
			sys.exit()
		if check_empty:check_file_empty(fl)

def check_file_empty(*files_in):
	"""
	Check if file is empty.

	Parameters
	----------
	files_in : str
		Input file (s) to check emptiness.

	Results
	-------
		return a boolean value.
	"""
	for fl in files_in:
		if not os.path.getsize(fl):
			logger.error(f"{fl} is empty.")
			sys.exit()

def get_file_prefix(file_in:str=None, times:int=1, split_symbol:str='.', include_path:bool=False):
	"""
	Get prefix of file, e.g. file is test.fa, then return test

	Parameters
	----------
	file_in : str
		Input file, relative path or abusolute path.
	times : int, default 1.
		How many times supposed to chop str
    split_symbol : str, default "."
        symbol to use to split file names
	include_path : boolean, default False
		Whether or not to include path.

	Returns
	-------
		Return prefix, include_path or not.
	"""
	#check_file_exist(file_in) #check outside, to escape file not exist but still want the prefix
	file_in = abs_path(file_in).split(split_symbol)
	file_in = split_symbol.join(file_in[0:len(file_in)-times])
	return include_path and file_in or os.path.basename(file_in)

def get_seqfile_prefix(seq_in:str=None):
	"""
	Get prefix of a sequence file.
	e.g. file is someting_1.fastq.gz, return something.

	Parameter
	---------
	seq_in : str
		Input file.

	Return
	------
	str:
		Return a string indicate as the prefix of seqfile.
	"""

	seq_in = os.path.basename(seq_in)
	return re.sub(".\d.fastq.gz|.\d.fq.gz|.\d.fastq|.\d.fq|.\d.fa.gz|.fasta.gz|.fa|.fasta", "", seq_in)

# Deprecated, to add times parameter
#	for i in range(times):
#		if include_path:return os.path.splitext(f)[0]
#		return os.path.splitext(os.path.basename(f))[0]

def get_file_path(*files_in):
	"""
	Get absolute path of input and return.

	Parameters
	----------
	files_in : str
		Input file (s)

	Return
	------
	str
		Absolute path (s) of your input file (s).
	"""
	final_paths = [os.path.dirname(abs_path(fl)) for fl in files_in]
	return len(final_paths) == 1 and final_path[0] or final_paths

def perfect_open(file_in:str=None):
	"""
	Make a perfect open for file

	Parameters
	----------
	filein : str
		input file to open

	Returns
	-------
	str
		Return file handle (s)
	"""
	return ".gz" in filein and gzip.open(filein, "rt") or open(filein, "r")

def close_file(*file_handles):
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
	for fh in file_handles:
		fh.close()

def find_files(dr:str=None, suffix='fa'):
	"""
	Find files with suffix.

	Parameters
	----------
	Dir : str
		Input directory to find files
	suffix : str, default "fa"
		suffix of files

	Returns
	-------
	list
		A list of found files with full path.
	"""
	return [dr+"/"+fl for fl in os.listdir(Dir) if fl.endswith(suffix)]

def parse_json(json_in):
	"""
	Parse json file.

	parameters
	----------
	json_in : str
		Input json format file.

	Returns
	-------
		Return a string but formatted like dataframe.
	"""
	import json
	with open(json_in) as fp:
	#	out.write('Head line.\n')
		json_file = json.load(fp)
	#	out_df = pd.DataFrame(index=None, columns=None)
		out_st_df = ""
		n0 = json_file["name"]
		for cld1 in json_file["children"]:
			n1 = cld1['name']
			for cld2 in cld1["children"]:
				n2 = cld2["name"]
				for cld3 in cld2["children"]:
					n3 = cld3["name"] ## "ko01000" is strange. what this for?
					if "children" not in cld3:
						#out.write('%s\t%s\t%s\t%s\t%s\n' %(n0, n1, n2, n3, n4))
						#out.write(f'{n0}\t{n1}\t{n2}\t{n3}\n')
						#out_df.append(pd.DataFrame([n0, n1, n2, n3]).T)
						out_st_df += f"{n0}\t{n1}\t{n2}\t{n3}\n"
					else:
	#					print(cld3)
						for cld4 in cld3["children"]:
							n4 = cld4["name"]
							if "children" not in cld4:
								out_st_df += f"{n0}\t{n1}\t{n2}\t{n3}\t{n4}\n"
							else:
								for cld5 in cld4["children"]:
									n5 = cld5["name"]
									if "children" not in cld5:
										out_st_df += f"{n0}\t{n1}\t{n2}\t{n3}\t{n4}\t{n5}\n"
									else:
										for cld6 in cld5["children"]
											n6 = cld6["name"]
											out_st_df += f"{n0}\t{n1}\t{n2}\t{n3}\t{n4}\t{n5}\t{n6}\n"
	return out_st_df
