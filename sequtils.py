
#####################################################
##												   ##
# seqs_utils.py -- fasta and fastq files utils	#
##												   ##
#####################################################

import Bio
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import gzip
import re

class gc:
	_gc_dict = {}
	def __init__(self, seqs, len_cutoff=0):
		self.seqs = seqs
		self.len_cutoff = len_cutoff
		#print(len_cutoff)
	def count_fasta_gc(self):
		with open(self.seqs) as fasta_handle:
			# use low-level parser to speed up when dealing with large data
			for title, seq in SimpleFastaParser(fasta_handle):
				if len(seq) >= self.len_cutoff:
					self._gc_dict[title] = [self._count_string_gc(seq), len(seq)]
		#	print(self._gc_dict)
			return self._gc_dict

	def count_fastq_gc(self):
		with open(self.seqs) as fastq_handle:
			# use low-level parser to speed up when dealing with large data
			for title, seq, qual in FastqGeneralIterator(fastq_handle):
				if len(seq) >= self.len_cutoff:
					self._gc_dict[title] = self._count_string_gc(seq)
			return self._gc_dict
	
	def _count_string_gc(self, string):
		_total_gc = string.count('G') + string.count('C')+string.count('g') + string.count('c')
		return round(_total_gc/len(string)*100., 2) # round to reduce float


def perfect_open(f):
	'''
	Make a perfect open for bio python, return a file handle.
	'''
#	if seqtype not in ['fastq', 'fastq']:
#		sys.exit('Wrong input file, must be fasta or fastq.')
	if re.search('gz$', f):
		seq_in = gzip.open(f, 'rt')
	else:
		seq_in = open(f, 'r')
	return seq_in

	


