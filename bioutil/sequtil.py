
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
from bioutil.system import files

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
					gc[t] = [cls._count_string_gc(seq)/len(seq)*100., len(seq)]
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
				gc[t] = cls._count_string_gc(seq)/len(seq)*100.
		fastq_handle.close()
		return gc
	
	def _count_string_gc(s):
		_total_gc = s.count('G') + s.count('C')+s.count('g') + s.count('c')
		return _total_gc

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

	@classmethod
	def evaluate_genome(cls, genome, ftype='file'):
		"""
		Evaluate genome and return genome features.

		Parameters:
		-----------
		genome:dict|file
			Input dict contains contigs or a FASTA file.
		ftype=dict|file
			input genome is dict or file [file]

		Returns:
		--------
		Return genome size, n50, maximal contig, minimal contig, gap number, gc ratio.
		"""
		if ftype != 'file' or ftype != 'dict':
			sys.exit('ftype can only be file or dict.')

		if ftype is file:
			genome = cls.read_fasta(genome)
		gap = 0
		gc = 0
		contig_lens = []
		for i in genome:
			contig = genome[i]
			contig_lens.append(len(contig))
			gap += len(re.findall('N+', contig))
			gc += cls._count_string_gc(contig) 

		genome_size = sum(contig_lens)
		gc = round(gc/genome_size*100., 2)

		contig_lens.sort(reverse=True)
		sum_len = 0
		for i in contig_lens:
			sum_len += i
			if sum_len >= genome_size*0.5:
				n50 = i
				break
		
		return genome_size, n50, max(contig_lens), min(contig_lens), gap, gc


class seqmodify:

	@classmethod
	def filter_fasta(cls, fasta, outfasta, 
				longer=None, shorter=None,
				first=0, end=0):
		"""
		Trim sequences according length.

		Parameters:
		-----------
		fasta:str
			input FASTA file.
		outfasta:str
			output FASTA file.
		longer=int
			exclude sequence longer than this cutoff, default None.
		shorter=int
			exclude sequence shorter than this cutoff, default None.
		first=int
			longest top n% sequences will be excluded, default 0.
		end=int
			shortest top n% sequences will be excluded, default 0.

		Returns:
		--------
		Trimmed FASTA
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
		-----------
		fasta:str
			fasta input file
		outf=int
			output length file, default is None.
		plot=bool
			plot a hist of fasta length, default is False. It has to set with outf
		bins=int
			number of bins to plot, default 10.
		
		Returns:
		--------
		Return dataframe contain FASTA length.
		"""
		
		length = {}
		in_handle = files.perfect_open(fasta)
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


	def extra_seq(fasta, idlist, outf=None, match=True):
		"""
		Extract FASTA sequence you need.
		
		Parameters:
		-----------
		fasta:str
			Input FASTA file.
		idlist:list or file contain a column of id, without header.
			idlist to extract corresponding sequences.
		outf=str
			output file name, default stout
		match=bool
			extract the sequence matched id or in turn.

		Result:
		--------
		Output FASTA file or stdout FASTA.

		"""
		if type(idlist) is str:
			idlist = pd.read_csv(idlist, sep='\t', header=None, index_col=0).index
		in_handle = files.perfect_open(fasta)
		outseq = []
		for rec in SeqIO.parse(in_handle, 'fasta'):
			if match:
				if str(rec.id) in idlist:
					outseq.append(rec)
			else:
				if str(rec.id) not in idlist:
					outseq.append(rec)
		if outf:
			SeqIO.write(outseq, outf, 'fasta')


	@classmethod
	def break_fasta(fasta, outfasta, symbol='N', exact=True):
		"""
		Use this function to break sequence using symbol (e.g. Ns).

		Parameters:
		-----------
		fasta:str
			Input FASTA file.
		outfasta:str
			Output FASTA file.
		symbol=str
			symbol to use to break FASTA. [N]
		exact=bool
			exact symbol or not, e.g, set symbol is NN, exact=True will not NNN or NNNN as a cut point. [True]

		Result:
		-------
		Output a FASTA file with broken using symbol.
		"""
	
		symbol_len = len(symbol)
		symbol += '+' # make a 're' match to indicate one or more symbol
		print(symbol)
		fasta = cls.read_fasta(fasta)
		c = 0
		start = 0
		end = 0
		with open(outfasta, 'w') as out:
			for i in fasta:
				gaps = re.findall(symbol, fasta[i])
				for gap in gaps:
					pos = fasta[i][end:].find(gap)
					end += pos
					## use symbol_len to replace, to judge whether to stop here or not.
					if len(gap) == symbol_len:
						c += 1
						out.write('>%s_%d|size=%s\n%s\n' % (i, c, end-start, fasta[i][start:end]))
						start = end + len(gap)
						end += len(gap)
						continue
					if exact:
						end += len(gap)
					else:
						c += 1
						out.write('>%s_%d|size=%s\n%s\n' % (i, c, end-start, fasta[i][start:end]))
						start = end + len(gap)
						end += len(gap)
				## make the last one, cause n gaps will chunk sequences into n+1 small sequences.
				out.write('>%s_%d|size=%s\n%s\n' % (i, c, len(fasta[i])-start, fasta[i][start:]))


# for test

if __name__== '__main__':
	import sys
	fasta = sys.argv[1]
#	seqmodify.break_fasta(fasta, aaaa)
	seqmodify.filter_fasta(aaaa, fasta+'.more500', shorter=500)

