
#####################################################
#												   	#
# seqsutils.py -- fasta and fastq files utils		#
#												   	#
#####################################################

from re import findall
import pandas as pd
from biosut.biosys import files

class sequtil:

	@classmethod
	def cal_seq_gc(cls, seq, length:bool=False):
		"""
		Count sequence gc ratio and length.
		
		Parameters:
		-----------
		seq:str
			FASTA/FASTQ(.gz) file
		length:bool
			return length or not. default False.
		len_cutoff=int
			sequence below this length wont include, default 0 means count all sequences.

		Returns:
		--------
			Return a dict with key-value pairs are seqid and list[gc ratio]
		"""
		gc = {}
		# use perfect_open to deal with*.gz files
		fh = files.perfect_open(seqs)
		# use low-level parser to speed up when dealing with super large data
		# Jie, 2020-06-23, use Heng Li's readfq instead, roughly, 15% slower than Bio,
		# it's acceptable while considering file size.
		for t, seq, _ in cls.seq_reader(fh):
			gc[t] = [cls._count_string_gc(seq)/len(seq)*100., len(seq)]
		fh.close()
		return gc
	
	def _count_string_gc(s):
		s = s.upper()
		return s.count('G') + s.count('C')

	# this is a copy-and-paste from https://github.com/lh3/readfq/blob/master/readfq.py
	def seq_reader(fh):
		"""
		sequence generator.
		fh:str
			file handle.
		"""
		last = None # this is a buffer keeping the last unprocessed line
		while True: # mimic closure; is it a bad idea?
			if not last: # the first record or a record following a fastq
				for l in fp: # search for the start of the next record
					if l[0] in '>@': # fasta/q header line
						last = l[:-1] # save this line
						break
			if not last: break
			#name, seqs, last = last[1:].partition(" ")[0], [], None
			name, seqs, last = last[1:], [], None # Jie, modified to keep the comment of id.
			for l in fp: # read the sequence
				if l[0] in '@+>':
					last = l[:-1]
					break
				seqs.append(l[:-1])
			if not last or last[0] != '+': # this is a fasta record
				yield name, ''.join(seqs), None  # yield a fasta record
				if not last: break
			else: # this is a fastq record
				seq, leng, seqs = ''.join(seqs), 0, []
				for l in fp: # read the quality
					seqs.append(l[:-1])
					leng += len(l) - 1
					if leng >= len(seq): # have read enough quality
						last = None
						yield name, seq, ''.join(seqs) # yield a fastq record
						break
				if last: # reach EOF before reading enough quality
					yield name, seq, None # yield a fasta record instead
					break
	
	@classmethod
	def read_seq(cls, fl : str, length : bool=False, qual : bool=False):
		"""
		Read fasta format file in.
		Parameters:
		-----------
		fl:str
			FASTA/FASTQ(.gz) file
		length:bool
			output length instead of sequence, default False.
		
		Returns:
		--------
			Return a dict contain seq id as key and seq as value.
		"""

		seqs = {}
		fh = files.perfect_open(fasta)
		for t, seq, _ in seq_reader(fh):
			seqs[t] = seq
			if length:seqs[t] = len(seq)
		fh.close()
		return seqs

	@classmethod
	def evaluate_genome(cls, genome):
		"""
		Evaluate genome and return genome features.

		Parameters
		----------
		genome : file
			Input file contains contigs or a FASTA file.

		Returns
		-------
			Return genome size, contig number, n50, maximal contig, minimal contig, gap number, gc ratio.
		"""
		fh = files.perfect_open(genome)
		gap, gc, contig_num, contig_len = 0, 0, 0, []
		for t, seq, _ in seq_reader(fh):
			contig_num += 1
			contig_len.append(len(seq))
			gap += len(findall('N+', seq))
			gc += cls._count_string_gc(seq) 
		
		genome_size = sum(contig_len)
		gc = round(gc/genome_size*100., 2)

		contig_len.sort(reverse=True)
		sum_len = 0
		for i in contig_len:
			sum_len += i
			if sum_len >= genome_size*0.5:
				n50 = i
				break
		
		return genome_size, contig_num, n50, max(contig_len), min(contig_len), gap, gc


class seqalter:
	@classmethod
	def filter_fasta(cls, fasta, outfasta, 
				longer : bool = None, shorter:bool = None,
				first = 0, end = 0):
		"""
		Trim sequences according length.

		Parameters
		----------
		fasta : str
			input FASTA file.
		outfasta : str
			output FASTA file.
		longer : int
			exclude sequence longer than this cutoff, default None.
		shorter : int
			exclude sequence shorter than this cutoff, default None.
		first : int or float
			longest top n% sequences will be excluded, default 0.
		end : int or float
			shortest top n% sequences will be excluded, default 0.

		Returns
		-------
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

	def cal_length(fasta, outf = None, plot : bool = False, bins = 10):
		"""
		calculate sequences length.
		
		Parameters
		----------
		fasta : str
			fasta input file
		outf : str, default None
			Output length file.
		plot : bool, default False
			Plot a hist of fasta length. It has to set with `outf`.
		bins : int, default is 10.
			Number of bins to plot.
		
		Returns
		-------
			Return dataframe contain FASTA length.
		"""
		length = {}
		fh = files.perfect_open(fasta)
		for t, seq in sequtil.seq_reader(fh):
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

	def extract_seq(in_file, idlist, in_type = 'fasta', outdir = None):
		"""
		Extract sequences you need.
		
		Parameters
		----------
		in_file : str
			Input sequence file.
		idlist : list, or file contain a column of id, file without header.
			idlist to extract corresponding sequences. id is the string before gap.
		in_type : str, default in_type=fasta
			input sequences type, fasta or fastq, default is fasta
		outdir:str
			output dir, default is the same as in_file directory.

		Result:
		-------
			Return matched idlist sequences and negtive matched idlist sequences.
		"""
		
		if in_type != 'fasta' and in_type != 'fastq':
			sys.exit('in_type can only be fasta or fastq.')

		if type(idlist) is str:
			#if low_level:
			df_reader = pd.read_csv(idlist, sep='\t', header=None, index_col=0, iterator=True)
			idlist = []
			while True:
				try:
					chunk = df_reader.get_chunk(500000)
					idlist.extend(list(chunk.index))
				except StopIteration:
					print("Finished looping the idlist in.")
					break
		
		if outdir:
			outdir = path.get_path(in_file)

		prefix = outdir + '/' + files.get_prefix(in_file)

		fh = files.perfect_open(in_file)
		match, negmatch = {}, {}
		
		match_out = open(prefix + '.match.' + in_type, 'w')
		negmatch_out = open(prefix + '.negmatch.' + in_type, 'w')

		if in_type == 'fasta':
			for t, seq, _ in readseq(fh):
				if t.partition(' ')[0] in idlist:
					match_out.write('>%s\n%s\n' % (t, seq))
				else:
					negmatch_out.write('>%s\n%s\n' % (t, seq))
		else:
			for t, seq, qual in readseq(fh):
				if t.partition(' ')[0] in idlist:
					match_out.write('@%s\n%s\n+\n%s\n' % (t, seq, qual))
				else:
					negmatch_out.write('@%s\n%s\n+\n%s\n' % (t, seq, qual))

		fh.close()
		match_out.close()
		negmatch_out.close()

	def split_fasta(fasta, outfasta, symbol = 'N', exact : bool = True):
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
		fh = files.perfect_open(fasta)
		print(symbol, symbol_len)
		with open(outfasta, 'w') as out:
			for t, seq, _ in sequtil.seq_reader(fh):
				c, start, end = 0, 0, 0
				gaps = findall(symbol, seq)

				if len(gaps) == 0:
					out.write('>%s\n%s\n' % (i, seq))
					continue

				for gap in gaps:
					pos = seq[end:].find(gap)
					end += pos
					## use symbol_len to replace, to judge whether to stop here or not.
					if len(gap) == symbol_len: # exact a 'gap' to split fasta
						c += 1
						out.write('>%s_%d|len=%s\n%s\n' % (i, c, end-start, seq[start:end]))
						start = end + len(gap)
						end += len(gap)
						continue
					if exact: # N is more than expected.
						end += len(gap)
					else:
						c += 1
						out.write('>%s_%d|len=%s\n%s\n' % (i, c, end-start, seq[start:end]))
						start = end + len(gap)
						end += len(gap)
				## make the last one, cause n gaps will chunk sequences into n+1 small sequences.
				out.write('>%s_%d|len=%s\n%s\n' % (i, c+1, len(seq)-start, seq[start:]))

# for test
if __name__== '__main__':
	import sys
	fasta = sys.argv[1]
#	seqalter.break_fasta(fasta, aaaa)
	seqalter.filter_fasta(aaaa, fasta+'.more500', shorter=500)

