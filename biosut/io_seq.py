"""
The :mod:`biosut.io_seq` includes utilities to operate sequence files.
"""

# Author: Jie Li <mm.jlli6t@gmail.com>
# License: GNU v3.0

import oss
from re import findall
from .gt_file import perfect_open

# copy-and-paste from https://github.com/lh3/readfq/blob/master/readfq.py
def iterator(fh, chop_comment:bool=False):
	"""
	Sequence iterator.
	fh : str
		Input file handle.
	chop_comment : bool, default is False
		Chop comment in sequence id or not.
	"""
	last = None # this is a buffer keeping the last unprocessed line
	while True: # mimic closure; is it a bad idea?
		if not last: # the first record or a record following a fastq
			for l in fh: # search for the start of the next record
				if l[0] in '>@': # fasta/q header line
					last = l[:-1] # save this line
					break
		if not last: break
		#name, seqs, last = last[1:], [], None # jlli6t, keep comment of seq id.
		name = last[1:].partition(" ")[0] if chop_comment else last[1:]
		seqs, last = [], None
		for l in fh: # read the sequence
			if l[0] in '@+>':
				last = l[:-1]
				break
			seqs.append(l[:-1])
		if not last or last[0] != '+': # this is a fasta record
			yield name, ''.join(seqs), None  # yield a fasta record
			if not last: break
		else: # this is a fastq record
			seq, leng, seqs = ''.join(seqs), 0, []
			for l in fh: # read the quality
				seqs.append(l[:-1])
				leng += len(l) - 1
				if leng >= len(seq): # have read enough quality
					last = None
					yield name, seq, ''.join(seqs) # yield a fastq record
					break
			if last: # reach EOF before reading enough quality
				yield name, seq, None # yield a fasta record instead
				break

def fq2fa(infq, outfa):
	"""
	Convert FASTQ format to FASTA format file.

	Parameters
	----------
	infq : str
		Input FASTQ(.gz) file.
	outfa : str
		Output FASTA file.

	Result
	------
		Output converted FASTA file.
	"""
	fh = perfect_open(fq)
	with open(outfa, 'w') as outf:
		for t, seq, _ in iterator(fh):
			outf.write('>%s\n%s\n' % (t, seq))
	fh.close()

def string_gc(string):
	"""
	Count string G/C number.

	Parameters
	----------
	string : str

	Returns
	-------
		Int :
            Return number of GC count.
	"""
	string = string.upper()
	return string.count('G') + string.count('C')

def gc_to_dict(inseq:str, len_cutoff:int=0, length:bool=False):
	"""
	Count GC  of sequences and other characteristics of sequences \
	to dict.

	Parameters
	----------
	inseq : str
		FASTA/FASTQ(.gz) file
	len_cutoff : int, default 0.
		Sequence below this length will be excluded.
	length : bool, default False.
		Also return length of each sequences or not.

	Returns
	-------
	dict :
        A dict with sequence id as key and along gc and \
		other characteristics of sequences as value.
	"""
	gc = {}
	# use perfect_open to deal with*.gz files
	fh = perfect_open(inseq)
	# use low-level parser to speed up when dealing with super large data
	# jlli6t, 2020-06-23, use Heng Li's readfq instead, roughly, 15% slower than Bio,
	# it's acceptable while considering file size.
	for t, seq, _ in iterator(fh):
		if len(seq)<len_cutoff:continue
		gc[t] = [string_gc(seq)]
		if length:gc[t].append(len(seq))
	fh.close()
	return gc

def seq_to_dict(inseq:str, outqual:bool=False, len_cutoff:int=0):
	"""
	Read and return sequences to dict format.

	Parameters
	----------
	inseq : str
		FASTA/FASTQ(.gz) file
	outqual : bool, default False.
		Include qual in output or not.
	len_cutoff : int
		Sequences below this cutoff will be discarded.
	Returns
	-------
		Return a dict contain seq id as key and seq as value.
	"""

	seqs = {}
	fh = perfect_open(inseq)
	for t, seq, _ in iterator(fh):
		if len(seq) < len_cutoff:continue
		seqs[t] = [seq]
		if outqual:seqs[t].append(_)
	fh.close()
	return seqs

def evaluate_genome(genome, len_cutoff:int=500):
	"""
	Evaluate genome and return genome traits.

	Parameters
	----------
	genome : file
		Input file contains contigs or a FASTA file.
	len_cutoff : int, default 500
		sequences below this length will be excluded.

	Returns
	-------
		Return genome size, contig number, n50, maximal contig, \
		minimal contig, gap number, gc ratio.
		"""

	fh = perfect_open(genome)
	gap, gc, contig_num, contig_len = 0, 0, 0, []
	for t, seq, _ in iterator(fh):
		if len(seq) < len_cutoff:continue
		contig_num += 1
		contig_len.append(len(seq))
		gap += len(findall('N+', seq))
		gc += string_gc(seq)

	genome_size = sum(contig_len)
	gc = round(gc/genome_size*100., 2)

	contig_len.sort(reverse=True)
	sum_len = 0
	for i in contig_len:
		sum_len += i
		if sum_len >= genome_size*0.5:
			n50 = i
			break
	return genome_size, contig_num, n50, max(contig_len), \
			min(contig_len), gap, gc
