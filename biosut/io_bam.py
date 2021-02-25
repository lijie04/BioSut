"""
The :mod:`biosut.io_bam` includes some bam operation.
"""

# Author: Jie Li <mm.jlli6t@gmail.com>
# License: GNU v3.0
# Copyrigth: 2015 -

import os
import pysam as ps

from .gt_exe import exe_cmd

# code relies on samtools, so please add samtools in.
def sort_bam(bam, overlay:bool=False):
	"""
	Sort bam file.

	Parameters
	----------
	bam : str
		Input bam file.
	overlay : bool, default False
		Overlay existing result or not.

	Returns
	-------
	str
		Sorted bam file.
	"""
	srt_bam = bam + '.srt.bam'
	if overlay:
		ps.sort('-o', srt_bam, bam)
		return srt_bam
	if is_sorted(bam):return bam
	ps.sort('-o', srt_bam, bam)
	return srt_bam

def index_bam(bam, overlay:bool=False):
	"""
	Index bam.

	Parameters
	----------
	bam : str
		Input bam file
	overlay : bool, default False
		overlay existing results or not.
	"""

	if overlay:
		if is_sorted(bam):ps.index(bam)
		return bam
		srt_bam = sort_bam(bam, overlay=True)
		ps.index(srt_bam)
		return srt_bam

	if os.path.isfile(bam+'.bai'):return bam

	try:
		ps.index(bam)
		return bam
	except ps.utils.SamtoolsError:
		srt_bam = sort_bam(bam)
		ps.index(srt_bam)
		return srt_bam

# deprecated use index_bam(overlay=True) instead
#	def sure_index_bam(bam):
#		"""
#		Indexing bam file no matter it is existed or not.

#		Parameters
#		----------
#		bam:str
#			input bam file.

#		Results
#		-------
#		Generate indexed bam file and return bam.
#		"""

#		srt_bam = bam + '.srt.bam'
#		try:
#			ps.index(bam)
#			return bam
#		except ps.utils.SamtoolsError: # bam is not sorted.
#			if not os.path.isfile(srt_bam):
#				cls.sort_bam(bam)
#			ps.index(srt_bam)
#			return srt_bam

def is_sorted(bam):
	"""
	Check a bam is sorted or not.

	Parameters
	----------
	bam : str
		input bam file.

	Returns
	-------
		Bool value.
	"""
	cmd = f'samtools view -H {bam}|head -n 1'
	out, err = exe_cmd(cmd)
	if 'unsorted' in out.decode():return False
	return True
