"""
The :mod:`biosut.alter_seq` includes utilities to operate sequence files.
"""

# Author: Jie Li <mm.jlli6t@gmail.com>
# License: GNU v3.0
# Copyrigth: 2015 -

from re import findall
import gzip

from . import io_seq
from .gt_file import perfect_open, close_file, get_file_prefix

def select_seq(inseq:str=None, outseq:str=None, longer:int=0, shorter:int=0, \
				first:float=0, end:float=0, outqual=False):
	"""
	Select sequences according length.

	Parameters
	----------
	infasta : str
		input FASTA/FASTQ(.gz) file.
	outseq : str
		output seq file in FASTA/FASTQ format.
	longer : int, default keep all sequence.
		exclude sequence longer than this cutoff.
	shorter : int
		exclude sequence shorter than this cutoff.
	first : float, 1 means all.
		longest top n% sequences will be excluded, default 0.
	end : float, 1 means all
		shortest top n% sequences will be excluded, default 0.
	outqual : bool, default False
		Include qual in output or not.

	Returns
	-------
		Trimmed FASTA/FASTQ sequences.
	"""
	fh = perfect_open(inseq)
	all_length = []
	for t, seq, _ in io_seq.iterator(fh):
		if shorter and len(seq) < shorter:continue
		if longer and len(seq) > longer:continue
		all_length.append(len(seq))
	fh.close()

	all_length.sort(reverse=True)
	total = len(all_length)
	if first:first = all_length[round(first * total + 0.5)]
	if end:end = all_length[total-round(end * total + 0.5)-1]
	fh = perfect_open(infasta)
	with open(outseq, 'w') as outf:
		for t, seq, qual in io_seq.iterator(fh):
			if first and len(seq) > first:continue
			if end and len(seq) < end:continue
			if outqual:
				outf.write(f'@{t}\n{seq}\n+\n{qual}\n')
			else:
				outf.write(f'>{t}\n{seq}\n')
	fh.close()

def split_fasta(infasta, outfasta, symbol:str='N', exact:bool=True):
	"""
	Split sequences using symbol (e.g. Ns).

	Parameters
	----------
	infasta : str
		Input FASTA file.
	outfasta : str
		Output FASTA file.
	symbol : str, default 'N'
		symbol to use to break sequence
	exact : bool, default True
		exact symbol or not, \
		e.g, set symbol to NN, and exact=True, \
		program will not recognize NNN or NNNN as a split site.

	Result
	------
		Output splitted sequence file.
	"""

	symbol_len = len(symbol)
	symbol += '+' # make a 're' match to indicate one or more symbol
	fh = perfect_open(infasta)
	out = open(outfasta, 'w')
	print(symbol, symbol_len)
	for t, seq, _ in io_seq.iterator(fh):
		c, start, end = 0, 0, 0
		gaps = findall(symbol, seq)

		if len(gaps) == 0:
			out.write('>%s\n%s\n' % (t, seq))
			continue

		for gap in gaps:
			pos = seq[end:].find(gap)
			end += pos
			# use symbol_len to replace, judge whether to stop here or not.
			if len(gap) == symbol_len: # exact a 'gap' to split fasta
				c += 1
				out.write(f'>{t}_{c}|len={end-start}\n{seq[start:end]}\n')
				start = end + len(gap)
				end += len(gap)
				continue
			if exact: # N is more than expected.
				end += len(gap)
				continue
			c += 1
			out.write(f'>{t}_{c}|len={end-start}\n{seq[start:end]}\n')
			start = end + len(gap)
			end += len(gap)
		# output the last one, as n gaps cut sequences into n+1 sequences.
		out.write(f'>{t}_{c+1}|len={len(seq)-start}\n{seq[start:]}\n')
	close_file(fh, out)

def reorder_PE_fq(infq1, infq2, outdir=None):
	"""
	Reorder pair-end FASTQ files to make fq1 & fq2 in same order.

	Parameters
	----------
	infq1 : str
		Input FASTQ 1 file (.gz)
	infq2 : str
		Input FASTQ 2 file (.gz)
	outdir : str, default None
		Outdir to output reordered files, without outdir, \
		old files will be replaced.

	Results
	-------
		Output pair-end FASTQ files that contain sequences in same order.
	"""
	fq1 = io_seq.seq_to_dict(infq1, qual=True, len_cutoff=0)
	fq2 = io_seq.seq_to_dict(infq2, qual=True, len_cutoff=0)

	if '.gz' in infq1:infq1 = get_file_prefix(infq1, include_path=True)
	if '.gz' in infq2:infq2 = get_file_prefix(infq2, include_path=True)

	if outdir:
		fq1_out = open(f'{outdir}/{os.path.basename(infq1)}', 'w')
		fq2_out = open(f'{outdir}/{os.path.basename(infq2)}', 'w')
	else:
		fq1_out = open(infq1, 'w')
		fq2_out = open(infq2, 'w')

	#	print('Reordering fastq sequences id.')
	for sid in fq1.keys():
		fq1_out.write(f'@\n{sid}\n{fq1[t][0]}\n+\n{fq1[t][1]}\n')
		fq2_out.write(f'@\n{sid}\n{fq2[t][0]}\n+\n{fq2[t][1]}\n')

	close_file(fq1, fq2)

def extract_seq(inseq:str=None, idlist:str=None, outseq:str=None, outqual:bool=False, \
				outseq2:str=None):
	"""
	Extract sequences in idlist.

	Parameters
	----------
	inseq : str
		Input FASTA/FASTQ(.gz) file.
	idlist : list
		idlist to extract corresponding sequences.
	outseq : str
		File to output matched sequences.
	outqual : boolean, default False
		Include qual in output or not if input sequence is in FASTQ format.
	outseq2 : str, default None
		If specified, negtive match sequences will be write into this file.

	Result
	------
		Output matched (and negtive matched) sequence file,\
		and return matched and negmatched sequences dicts.
	"""
#		all_id = []
#		if type(idlist) is str:
#			with open(idlist) as id_in:
#				for Id in id_in:
#					all_id.append(Id.strip())
#		else:
#			all_id = idlist

	match, negmatch = {}, {}
	match_out = open(outseq, 'w')
	if out_negmatch:negmatch_out = open(outseq + '.negmatch', 'w')
	fh = perfect_open(inseq)

	if outqual:
		for t, seq, _ in io_seq.iterator(fh):
			if t in idlist:
				match_out.write(f'@{t}\n{seq}\n+\n{_}\n')
				match[t] = [seq, _]
				continue
			negmatch[t] = [seq, _]
			if out_negmatch:negmatch_out.write(f'@{t}\n{seq}\n+\n{_}\n')
	else:
		for t, seq, _ in io_seq.iterator(fh):
			if t in idlist:
				match_out.write(f'>{t}\n{seq}\n')
				match[t] = [seq]
				continue
			negmatch[t] = [seq]
			if out_negmatch:negmatch_out.write(f'>{t}\n{seq}\n')
	close_file(fh, match_out)
	if out_negmatch:negmatch_out.close()
	return match, negmatch

def trim_headn(inseq:str=None, outseq:str=None, outqual:bool=False):
	"""
	Trim N from head of seq.

	parameter
	---------
	inseq : str
		input sequence FASTA/FASTQ(.gz) file.
	outseq : str
		output file name, endswith gz will compress output file.
	outqual : boolean, default is False.
		output quality or not if input sequence is in FASTQ format.

	Result
	------
	Output trimmed N sequences into output file.
	"""

	def remove_first_N(string):
		while string[0] in 'Nn':
			string = string[1:]
		return string

	if not inseq:
		logger.info('Have to specify an seq file for trim_headn.')
		sys.exit()
	if not outseq:
		logger.info('Have to specify an seq file for trim_headn.')
		sys.exit()

	fh = perfect_open(inseq)
	oh = '.gz' in outseq and gzip.open(outseq, 'wb') or open(outseq, 'w')
	for t, seq, _ in io_seq.iterator(fh):
		seq = remove_first_N(seq)
		line = outqual and f'@{t}\n{seq}\n+\n{_[len(_)-len(seq):]}\n' or f'>{t}\n{seq}\n'
		oh.write('.gz' in outseq and line.encode() or line)
	close_file(fh, oh)
