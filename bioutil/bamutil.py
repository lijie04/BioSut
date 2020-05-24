
#####################################################################
##																	#
# bam_utils.py -- some bam file related util functions				#
##																	#
####################################################################


import os
from bioutil.system import path
from bioutil.sequtil import sequtil
import pysam as ps
from re import findall
import pandas as pd

class bamutil:

	@classmethod
	def index_bam(cls, bam, override=False):
		"""
		Indexing bam file.

		Parameters:
		bam:str
			input bam file
		override:bool
			True/False to override existing results, default True.
		"""

		if override:
			bam = cls._sure_index_bam(bam)
			return bam

		if not os.path.isfile(bam + '.bai'):
			try:
				ps.index(bam)
				return bam
			except ps.utils.SamtoolsError:
				if not os.path.isfile(srt_bam):
					srt_bam = cls.sort_bam(bam, override=override)
					ps.index(srt_bam)
					return srt_bam
				else:
					ps.index(srt_bam)
					return srt_bam
		else:
			return bam
	
	def _sure_index_bam(bam):
		"""
		Indexing bam file no matter it is existed or not.

		Parameters:
		-----------
		bam:str
			input bam file.

		Results:
		--------
		Generate indexed bam file and return bam.
		"""
		
		srt_bam = bam + '.srt.bam'
		try:
			ps.index(bam)
			return bam
		except ps.utils.SamtoolsError: # indicate bam is not sorted.
			if not os.path.isfile(srt_bam):
				cls.sort_bam(bam)
			ps.index(srt_bam)
			return srt_bam


	def sort_bam(bam, override=False):
		"""
		Sort bam file.

		Parameters:
		bam:str
			input bam file.
		override:bool
			True/False to override existing result, default is True.

		Returns:
		str
			final bam file name
		"""

		srt_bam = bam + '.srt.bam'
		if override:
			ps.sort('-o', srt_bam, bam)
			return srt_bam
		# don't know how to judge a bam is sorted or not, so judge with name which is my naming style.
		if not os.path.isfile(srt_bam):
			ps.sort('-o', srt_bam, bam)
		return srt_bam
	
	@classmethod
	def retrieve_reads(cls, bam, ref, out_prefix, **kargs):
		"""
		Extract reads from bam according the reference provided.

		Parameters:
		bam:str
			bam file to process
		ref:str
			referece fasta file to use
		out_prefix:str
			output prefix to use for output files, must include paths.
		
		secondary=bool
			keep secondary alignments.
		qcfail=bool
			set to keep qcfail alignments according to samtools principle.
		duplicates=bool
			set to keep duplicates alignments.
		supplementary=bool
			set to keep supplementary alignments.
		reorder=bool
			set to reorder paired sequences, as they are not in same order originally.
		
		Results:
			Out put extracted reads into paired 1&2 and unpaired 3 types.
			out_prefix.1.fastq.gz, out_prefix.2.fastq.gz, out_prefix.forward.fastq.gz, out_prefix.reverse.fastq.gz and stats.xls file to tell you how many reads you have
			## And return paired fastq 1 &2, unpaired forward fastq 1 and reverse fastq 2 file.
			## Not return fq files at this moment.
		"""

		flag = {
				'secondary': '-F 0X100 ',
				'qcfail': '-F 0X200 ',
				'duplicates': '-F 0X400 ',
				'supplementary': '-F 0X800 '
				}

		cls.index_bam(bam)
		ref = ' '.join(sequtil.read_fasta(ref).keys())
		fq1 = '%s.1.fastq' % out_prefix
		fq2 = '%s.2.fastq' % out_prefix
		fq_forward = '%s.forward.fastq' % out_prefix
		fq_reverse = '%s.reverse.fastq' % out_prefix

		samtools_info = []

		cmd = 'samtools view -b '
		if kargs['secondary']:cmd += flag['secondary']
		if kargs['qcfail']:cmd += flag['qcfail']
		if kargs['duplicates']: cmd += flag['duplicates']
		if kargs['supplementary']:cmd += flag['supplementary']
		cmd_pair = cmd + '-f 0X2 %s %s|samtools fastq -1 %s -2 %s' % \
					(bam, ref, fq1, fq2)
		out, err = path.exe_proc(cmd_pair)
		samtools_info.append(err)
		cls._re_order_pair_fq(fq1, fq2)


		cmd_forward = cmd + '-F 0X2 -F OX10 %s %s|samtools fastq -o %s -0 %s' % \
						(bam, ref, fq_forward, fq_forward)
		out, err = path.exe_proc(cmd_forward)
		samtools_info.append(err)

		cmd_reverse = cmd + '-F 0X2 -f 0X10 %s %s|samtools fastq -o %s -0 %s' % \
						(bam, ref, fq_reverse, fq_reverse)
		out, err = path.exe_proc(cmd_reverse)
		samtools_info.append(err)

		cmd = 'gzip -f %s -2 %s %s %s' % \
				(fq1, fq2, fq_forward, fq_reverse)
		path.exe_proc(cmd)
		stats = cls._parse_samtools_info(samtools_info)
		stats.to_csv(out_prefix + '.retreived_reads.stat.xls', sep='\t')
		return stats

	def _parse_samtools_info(errs):
		stats = {}
		stats['singletons'] = {}
		stats['retreived reads'] = {}
		flags = ['paired', 'forward', 'reverse']
		
		for flg, err in zip(flags, errs):
			err = findall('discarded (\d+) singletons\\n.*processed (\d+) reads', err.decode())[0]
			stats['singletons'][flg] = err[0]
			stats['retreived reads'][flg] = err[1]
		
		return pd.DataFrame.from_dict(stats).astype(int)


	def _re_order_pair_fq(fq1, fq2):
		fq1_out = fq1
		fq2_out = fq2
		fq1 = sequtil.read_fastq(fq1, qual=True)
		fq2 = sequtil.read_fastq(fq2, qual=True)
	#	print('Reordering fastq sequences id.')
		with open(fq1_out, 'w') as fq1_out, open(fq2_out, 'w') as fq2_out:
			for seq_id in fq1.keys():
				fq1_out.write('@' + seq_id + '\n' + fq1[seq_id][0] + '\n+\n' + fq1[seq_id][1] + '\n')
				fq2_out.write('@' + seq_id + '\n' + fq2[seq_id][0] + '\n+\n' + fq2[seq_id][1] + '\n')


def infer_identity(aln):
	""" Infer identity of alignment, denominator not include hard-clipped bases, soft-clipped bases are included, I consider soft-clipp as gap-open"""
	matched_bases = [int(i) for i in findall(r"(\d*)M", aln.cigarstring)] ## M indicates matched, so extracted all "Matched" bases
	return sum(matched_bases)*100 / aln.infer_query_length()  ## use infer_query_length() as query may been trimmed while mapping (not trimming while QC)


def infer_alignlen_ratio(aln):
	"""Infer alignment length ratio weighting by Denominator has hard-clipped bases included to infer the real read length"""
	return aln.query_alignment_length*100 / aln.infer_read_length() ## here I have hard-clipped bases included

#def judge_insert(aln):

#def judge_edge_alignment(aln):
