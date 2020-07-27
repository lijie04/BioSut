"""
The :mod:`biosut.bioseq` includes utilities to operate sequence files.
"""
# Author: Jie Li <mm.jlli6t@gmail.com>
# License: GNU v3.0

from re import findall
from .biosys import gt_file

class io_seq:

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

	@classmethod
	def fq2fa(cls, infq, outfa):
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

		fh = gt_file.perfect_open(fq)
		with open(outfa, 'w') as outf:
			for t, seq, _ in cls.seq_reader(fh):
				outf.write('>' + t + '\n' + seq + '\n')
		fh.close()

	def string_gc(string):
		"""
		Count string G/C number.

		Parameters
		----------
		string : str

		Returns
		-------
			Return number of GC count.

		Examples
		--------
		>>> from biosut.bioseq.io_seq iport string_gc
		>>> string = "GATDGGBKEREWCGSGCGEW"
		>>> gc = string_gc(string)
		>>> gc
		7
		"""
		string = string.upper()
		return string.count('G') + string.count('C')

	@classmethod
	def gc_to_dict(cls, inseq:str, len_cutoff:int = 0, length:bool = False):
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
			Return a dict with sequence id as key and along gc and \
			other characteristics of sequences as value.
		"""
		gc = {}
		# use perfect_open to deal with*.gz files
		fh = gt_file.perfect_open(seqs)
		# use low-level parser to speed up when dealing with super large data
		# jlli6t, 2020-06-23, use Heng Li's readfq instead, roughly, 15% slower than Bio,
		# it's acceptable while considering file size.
		for t, seq, _ in cls.iterator(fh):
			if len(seq)<len_cutoff:continue
			gc[t] = [cls.string_gc(seq)]
			if length:gc[t].append(len(seq))
		fh.close()
		return gc

	@classmethod
	def seq_to_dict(cls, inseq:str, outqual:bool = False, len_cutoff:int=0):
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
		fh = gt_file.perfect_open(fasta)
		for t, seq, _ in iterator(fh):
			if len(seq) < len_cutoff:continue
			seqs[t] = [seq]
			if outqual:seqs[t].append(_)
		fh.close()
		return seqs

	@classmethod
	def evaluate_genome(cls, genome, len_cutoff:int=500):
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

		fh = gt_file.perfect_open(genome)
		gap, gc, contig_num, contig_len = 0, 0, 0, []
		for t, seq, _ in cls.iterator(fh):
			if len(seq) < len_cutoff:continue
			contig_num += 1
			contig_len.append(len(seq))
			gap += len(findall('N+', seq))
			gc += cls.string_gc(seq)

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

class alter_seq:
	@classmethod
	def select_seq(cls, inseq, outseq, longer = None, shorter = None, \
					first:float = 0, end:float = 0, outqual=False):
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
		fh = gt_file.perfect_open(infasta)
		all_length = []
		for t, seq, _ in io_seq.iterator(fh):
			if shorter and len(seq) < shorter:continue
			if longer and len(seq) > longer:continue
			all_length.append(len(seq))
		fh.close()

		all_length.sort(reverse=True)
		total = len(all_length)
		if first:first = all_length[round(first * total) + 0.5)]
		if end:end = all_length[total-round(end * total + 0.5)-1]
		fh = gt_file.perfect_open(infasta)
		with open(outseq, 'w') as outf:
			for t, seq, _ in io_seq.iterator(fh):
				if first and len(seq) > first:continue
				if end and len(seq) < end:continue
				if outqual:
					outf.write('@%s\n%s\n+\n%s\n'%(t, seq, _))
				else:
					outf.write('>%s\n%s\n'%(t, seq))
		fh.close()

	def split_fasta(infasta, outfasta, symbol = 'N', exact:bool = True):
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
		fh = gt_file.perfect_open(infasta)
		out = open(outfasta, 'w')
		print(symbol, symbol_len)
		for t, seq, _ in sequtil.seq_reader(fh):
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
					out.write('>%s_%d|len=%s\n%s\n' % \
								(t, c, end-start, seq[start:end]))
					start = end + len(gap)
					end += len(gap)
					continue
				if exact: # N is more than expected.
					end += len(gap)
					continue
				c += 1
				out.write('>%s_%d|len=%s\n%s\n' % \
				 			(t, c, end-start, seq[start:end]))
				start = end + len(gap)
				end += len(gap)
			# output the last one, as n gaps cut sequences into n+1 sequences.
			out.write('>%s_%d|len=%s\n%s\n' % \
						(t, c+1, len(seq)-start, seq[start:]))
		fh.close()
		out.close()

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

		if '.gz' in infq1:infq1 = gt_file.get_prefix(infq1, include_path=True)
		if '.gz' in infq2:infq2 = gt_file.get_prefix(infq2, include_path=True)

		if outdir:
			fq1_out = open('%s/%s' % (outdir, os.path.basename(infq1)), 'w')
			fq2_out = open('%s/%s' % (outdir, os.path.basename(infq2)), 'w')
		else:
			fq1_out = open(infq1, 'w')
			fq2_out = open(infq2, 'w')

	#	print('Reordering fastq sequences id.')
		for t in fq1.keys():
			fq1_out.write('@\n%s\n%s\n+\n%s\n' % (sid, fq1[t][0], fq1[t][1]))
			fq2_out.write('@\n%s\n%s\n+\n%s\n' % (sid, fq2[t][0], fq2[t][1]))

		fq1_out.close()
		fq2_out.close()

	def extract_seq(inseq, idlist, outseq, outqual:bool=False, \
					out_negmatch:bool=False):
		"""
		Extract sequences you need.

		Parameters
		----------
		inseq : str
			Input FASTA/FASTQ(.gz) file.
		idlist : list
			idlist to extract corresponding sequences.
		outseq : str
			File to output matched sequences.
		outqual : bool, default False
			Include qual in output or not
		out_negmatch : bool, default False
			Output negtive match sequences into *.negmatch or not.

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
		fh = gt_file.perfect_open(inseq)

		if outqual:
			for t, seq, _ in io_seq.iterator(fh):
				if t in idlist:
					match_out.write('@%s\n%s\n+\n%s\n' % (t, seq, _))
					match[t] = [seq, _]
					continue
				negmatch[t] = [seq, _]
				if out_negmatch:negmatch_out.write('@%s\n%s\n+\n%s\n' % \
				 									(t, seq, _))
		else:
			for t, seq, _ in io_seq.iterator(fh):
				if t in idlist:
					match_out.write('>%s\n%s\n' % (t, seq))
					match[t] = [seq]
					continue
				negmatch[t] = [seq]
				if out_negmatch:negmatch_out.write('>%s\n%s\n' % (t, seq))
		fh.close()
		match_out.close()
		negmatch_out.close()
		return match, negmatch
