#####################################################################
##																	#
# bam_utils.py -- some bam file related util functions				#
##																	#
####################################################################

import os
import pysam as ps
import subprocess as subp
from re import findall

class check:

	def check_samtools():
		"""Check if samtool is in your system path. If not existed, program exits. """
		# assum that a successful samtools returns, otherwise returns something non-zero
		try:
			subp.run(['samtools'], stdout=open(os.devnull, "w"), stderr=subp.STDOUT)
		except:
			logger.error(" [Error] Make sure samtools is in your system path.")
			sys.exit(" [Error] Seems samtools is not in your system path.")

	def check_bam_index(bam):
		"""First check if bam index exists, if not, then try to judge bam sorted or not,
	   	finally index bam and RETURN final bam file """

		if not os.path.isfile(bam +".bai"):
			try:
				ps.index(bam)
			except ps.utils.SamtoolsError:
				bam_prefix = os.path.splitext(bam)[0]
				if os.path.isfile(bam_prefix + ".srt.bam"):  # don't know how to judge a bam is sorted or not, so judge with name (which is my naming style)
					bam = bam_prefix + ".srt.bam"
					if not os.path.isfile(bam + ".bai"):
						ps.index(bam) # if find a bam that considered as sorted, then, index bam directly.
				else:
					#logger.warning(str(dt.datetime.now()) + " [WARNNING] Seems like your bam file is unsorted.")
					#logger.warning(str(dt.datetime.now()) + " [WARNNING] Therefore, I started to sort bam file.")
					ps.sort("-o", bam_prefix + ".srt.bam", bam)
					#logger.info(str(dt.datetime.now()) + " [INFO] I started to index bam.")
					bam = bam_prefix + "_srt.bam"
					ps.index(bam)
		return bam

def infer_identity(aln):
	""" Infer identity of alignment, denominator not include hard-clipped bases, soft-clipped bases are included, I consider soft-clipp as gap-open"""
	matched_bases = [int(i) for i in findall(r"(\d*)M", aln.cigarstring)] ## M indicates matched, so extracted all "Matched" bases
	return sum(matched_bases)*100 / aln.infer_query_length()  ## use infer_query_length() as query may been trimmed while mapping (not trimming while QC)


def infer_alignlen_ratio(aln):
	"""Infer alignment length ratio weighting by Denominator has hard-clipped bases included to infer the real read length"""
	return aln.query_alignment_length*100 / aln.infer_read_length() ## here I have hard-clipped bases included

def judge_insert(aln):

def judge_edge_alignment(aln):
