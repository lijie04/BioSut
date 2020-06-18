
#####################################################
#												   	#
# seqsutils.py -- fasta and fastq files utils		#
#												   	#
#####################################################

from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from re import findall
import pandas as pd
from bioutil.system import files
from bioutil.readseq import readseq

class sequtil:
	
	@classmethod
	def cal_fasta_gc(cls, seqs, len_cutoff=0):
		"""
		Count fasta file's sequence gc content and length.
		
		Parameters:
		seqs:	fasta file
		len_cutoff=cutoff of length to not count gc for this sequence. default 0.

		Result:
		Return a dict with key-value pairs are seqid and list[gc ratio, length of seq]
		"""
		gc = {}
		fh = files.perfect_open(seqs)
			# use low-level parser to speed up when dealing with super large data
		for t, seq in SimpleFastaParser(fh):
			if len(seq) >= len_cutoff:
				gc[t] = [cls._count_string_gc(seq)/len(seq)*100., len(seq)]
		fh.close()
		return gc
	
	@classmethod
	def cal_fastq_gc(cls, seqs, len_cutoff=0):
		"""
		Count fastq file's sequence gc content and length.
		
		Parameters:
		seqs:	fastq file
		len_cutoff=cutoff of length to not count gc for this sequence, default is 0, means will count all sequences' length.

		Result:
		Return a dict with key-value pairs are seqid and list[gc ratio, length of seq]
		"""

		gc = {}
		fh = files.perfect_open(seqs)
		# with open(seqs) as fastq_handle: this can't deal with *.gz file.
		# use low-level parser to speed up when dealing with large data
		for t, seq, _ in FastqGeneralIterator(fh):
			if len(seq) >= len_cutoff:
				gc[t] = cls._count_string_gc(seq)/len(seq)*100.
		fh.close()
		return gc
	
	def _count_string_gc(s):
		s = s.upper()
		return s.count('G') + s.count('C')

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
		fh = files.perfect_open(fasta)
		for t, seq in SimpleFastaParser(fh):
			seqs[t] = seq
			if length:seqs[t] = len(seq)
		fh.close()
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
		fh = files.perfect_open(fastq)
		# use lh3's code, more efficiently.
		# use low-level parser to deal with super large data, faster than lh3 readfq
		for t, seq, _ in FastqGeneralIterator(fh):
			seqs[t] = seq
			if length:seqs[t] = len(seq)
			if qual:seqs[t] = [seq, _]
		fh.close()
		return seqs

	@classmethod
	def evaluate_genome(cls, genome):
		"""
		Evaluate genome and return genome features.

		Parameters:
		-----------
		genome:file
			Input file contains contigs or a FASTA file.

		Returns:
		--------
		Return genome size, n50, maximal contig, minimal contig, gap number, gc ratio.
		"""
		fh = files.perfect_open(genome)
		gap, gc, contig_lens = 0, 0, []
		for t, seq in SimpleFastaParser(fh):
			contig_lens.append(len(seq))
			gap += len(findall('N+', seq))
			gc += cls._count_string_gc(seq) 

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

class seqalter:

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

		matched, _ = cls.extract_seq(fasta, length.index, in_type='fasta')
		SeqIO.write(matched, outfasta, 'fasta')


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
		fh = files.perfect_open(fasta)
		for t, seq in SimpleFastaParser(fh):
			length[t] = len(seq)
		new = {}
		new['length'] = length
		new = pd.DataFrame.from_dict(new)
		if outf:
			new.to_csv(outf, sep='\t')
			if plot:
				new.hist(column='length', grid=False, bins=bins)
				plt.savefig(outf+'.hist.pdf', dpi=600)
		return new

	def extract_seq(in_file, idlist, in_type='fasta',
					match_out=None, negmatch_out=None):
		"""
		Extract sequences you need.
		
		Parameters:
		-----------
		in_file:str
			Input sequence file.
		idlist:list or file contain a column of id, file without header.
			idlist to extract corresponding sequences. id is the string before gap.
		in_type=str
			input sequences type, fasta or fastq, default is fasta
		low_level:bool
			bool value to process as low level or not. default False
		match_out:str
			file to output positive matched items, default None.
		negmatch_out:str
			file to output negtive matched items, default None.

		Result:
		-------
			Return matched idlist sequences and negtive matched idlist sequences.
		"""

		if type(idlist) is str:
			#if low_level:
			df_reader = pd.read_csv(idlist, sep='\t', header=None, index_col=0, iterator=True)
			idlist = []
			while True:
				try:
					chunk = df_reader.get_chunk(10000000)
					idlist.extend(list(chunk.index))
				except StopIteration:
					print("Finished looping the db info in.")
					break
		if in_type != 'fasta' and in_type != 'fastq':
			sys.exit('in_type can only be fasta or fastq.')

		fh = files.perfect_open(in_file)
		match, negmatch = {}, {}
		if in_type == 'fasta':
			for t, seq in SimpleFastaParser(fh):
				if t.partition(' ')[0] in idlist:
					match[t] = seq
				else:
					negmatch[i] = seq
		else:
			for t, seq, qual in FastqGeneralIterator(fh):
				if t.partition(' ')[0] in idlist:
					match[t] = [seq, qual]
				else:
					negmatch[t] = [seq, qual]
		fh.close()

		if match_out:
			with open(match_out, 'w') as m_out:
				if in_type == 'fasta':
					for t in match:
						m_out.write('>', t, '\n', match[t], '\n')
				else:
					for t in match:
						m_out.write('@', t, '\n', match[t][0], '\n', '+', '\n', match[t][1], '\n')
		if negmatch_out:
			with open(negmatch_out, 'w') as nm_out:
				if in_type == 'fasta':
					for t in negmatch:
						nm_out.write('>', t, '\n', negmatch[t], '\n')
				else:
					for t in negmatch:
						nm_out.write('@', t, '\n', negmatch[t][0], '\n', '+', '\n', negmatch[t][1], '\n')

		return match, negmatch

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
		fasta = sequtil.read_fasta(fasta)
		print(symbol, symbol_len)
		with open(outfasta, 'w') as out:
			for i in fasta:
				c, start, end = 0, 0, 0
				contig = fasta[i]
				gaps = findall(symbol, contig)

				if len(gaps) == 0:
					out.write('>%s\n%s\n' % (i, contig))
					continue

				for gap in gaps:
					pos = contig[end:].find(gap)
					end += pos
					## use symbol_len to replace, to judge whether to stop here or not.
					if len(gap) == symbol_len:
						c += 1
						out.write('>%s_%d|len=%s\n%s\n' % (i, c, end-start, contig[start:end]))
						start = end + len(gap)
						end += len(gap)
						continue
					if exact:
						end += len(gap)
					else:
						c += 1
						out.write('>%s_%d|len=%s\n%s\n' % (i, c, end-start, contig[start:end]))
						start = end + len(gap)
						end += len(gap)
				## make the last one, cause n gaps will chunk sequences into n+1 small sequences.
				out.write('>%s_%d|len=%s\n%s\n' % (i, c+1, len(contig)-start, contig[start:]))


# for test

if __name__== '__main__':
	import sys
	fasta = sys.argv[1]
#	seqalter.break_fasta(fasta, aaaa)
	seqalter.filter_fasta(aaaa, fasta+'.more500', shorter=500)

