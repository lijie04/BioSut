
#####################################################
##												   ##
# seqsutils.py -- fasta and fastq files utils		#
##												   ##
#####################################################

import Bio
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import gzip
import re
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from system import files

class sequtil:
	
	@classmethod
	def count_fasta_gc(cls, seqs, len_cutoff=0):
		"""
		Count fasta file's sequence gc content and length.
		
		Parameters:
		seqs:	fasta file
		len_cutoff=cutoff of length to not count gc for this sequence. default 0.

		Result:
		Return a dict with key-value pairs are seqid and list[gc ratio, length of seq]
		"""
		gc = {}
		with open(self.seqs) as fasta_handle:
			# use low-level parser to speed up when dealing with super large data
			for t, seq in SimpleFastaParser(fasta_handle):
				if len(seq) >= len_cutoff:
					gc[t] = [cls._count_string_gc(seq), len(seq)]
		return gc
	
	@classmethod
	def count_fastq_gc(cls, seqs, len_cutoff=0):
		"""
		Count fastq file's sequence gc content and length.
		
		Parameters:
		seqs:	fastq file
		len_cutoff=cutoff of length to not count gc for this sequence, default is 0, means will count all sequences' length.

		Result:
		Return a dict with key-value pairs are seqid and list[gc ratio, length of seq]
		"""
		gc = {}
		fastq_handle = files.perfect_open(seqs)
		# with open(seqs) as fastq_handle: this can't deal with *.gz file.
		# use low-level parser to speed up when dealing with large data
		for t, seq, _ in FastqGeneralIterator(fastq_handle):
			if len(seq) >= len_cutoff:
				gc[t] = cls._count_string_gc(seq)
		fastq_handle.close()
		return gc
	
	def _count_string_gc(s):
		_total_gc = s.count('G') + s.count('C')+s.count('g') + s.count('c')
		return round(_total_gc/len(s)*100., 2) # round to reduce float

	def read_fasta(fasta, length=False):
		"""
		Read fasta format file in.
		Parameters:
		-----------
		fasta:str
			fasta format file
		length:bool
			output length instead of sequence, default False.
		
		Returns:
		--------
		Return a dict as id & sequence/length of seqeunce as key-value pairs.
		"""
		seqs = {}
		for rec in SeqIO.parse(fasta, 'fasta'):
			seqs[str(rec.id)] = str(rec.seq)
			if length:seqs[str(rec.id)] = len(str(rec.seq))
		return seqs

	def read_fastq(fastq, length=False, qual=False):
		"""
		Read fastq format file in.
		Parameters:
		-----------
		fastq:str
			fastq format file in
		length:bool
			output length.
		"""
		if length and qual:
			sys.exit('Cant obtain qual is along with sequences, not sequence length.')
		seqs = {}
		fastq_handle = files.perfect_open(fastq)
		# use low-level parser to deal with super large data
		for t, seq, _ in FastqGeneralIterator(fastq_handle):
			seqs[t] = seq
			if length:seqs[t] = len(seq)
			if qual:seqs[t] = [seq, _]
		fastq_handle.close()
		return seqs


class seqmodify:

	@classmethod
	def trimfa(cls, fasta, outfasta, 
				longer=None, shorter=None,
				first=0, end=0):
		"""
		Trim sequences according length.

		Parameters:
		-----------
		fasta:	input fasta file.
		outfasta:	output fasta file.
		longer=remove sequence longer than this cutoff.
		shorter=remove sequence shorter than this cutoff
		first=longest top n% sequences will be discard.
		end=shortest top n% sequences will be discard.

		Returns
		-------
		str
			Trimmed fa
		"""
		
		length = cls.cal_length(fasta)

		if longer:length = length.loc[length[length.length<=longer].index]
		if shorter:length = length.loc[length[length.length>=shorter].index]

		length = length.sort_values(by='length', ascending=False)
		first = round(first * len(length)/float(100) + 0.5)
		end = round(end * len(length)/float(100) + 0.5)
		length = length.iloc[first:len(length)-end, ]

		cls.extra_seq(fasta, length.index, outf=outfasta)

	def cal_length(fasta, outf=None, plot=False, bins=10):
		"""
		calculate sequences length.
		
		Parameters:
		fasta:	fasta input file
		outf=output length file, default is None.
		plot=bool, plot a hist of fasta length, default is False. It has to set with outf
		bins=number of bins to plot, default 10.
		"""
		
		length = {}
		in_handle = perfect_open(fasta)
		for rec in SeqIO.parse(in_handle, 'fasta'):
			length[str(rec.id)] = len(rec.seq)
		new = {}
		new['length'] = length
		new = pd.DataFrame.from_dict(new)
		if outf:
			new.to_csv(outf, sep='\t')
			if plot:
				new.hist(column='length', grid=False, bins=bins)
				plt.savefig(outf+'.hist.pdf', dpi=600)
		return new

	def extra_seq(fasta, idlist, outf=None):
		if type(idlist) is str:
			idlist = pd.read_csv(idlist, sep='\t', header=None, index_col=0).index
		in_handle = perfect_open(fasta)
		outseq = []
		for rec in SeqIO.parse(in_handle, 'fasta'):
			if str(rec.id) in idlist:
				outseq.append(rec)
		if outf:
			print('printing.')
			SeqIO.write(outseq, outf, 'fasta')


