"""
The :mod:`biosut._diamond` integrated diamond related operations.
"""

__author__ = 'Jie Li'
__copyright__ = 'Copyright 2023'
__credits__ = 'Jie Li'
__license__ = 'GPLv3.0'
__maintainer__ = 'Jie Li'
__email__ = 'jlli6t near gmail.com'

import os
from os.path import join
import pandas as pd
from numpy import unique
from loguru import logger

from biosut import biosys


class Aligner:
	def __init__(self, query, subject, outdir):
		self.query = query
		self.subject = subject
		self.outdir = outdir

	def align(self, query, d_type, outdir, arg):

		taxid = {'Archaea':'2157', 'Bacteria':'2'} # define tax id
		if not os.path.isfile(db + '.dmnd'): # check existance of dimaond index, if not, create index
			cls._index_diamond_db(db)

		biosys.check_file_exist(query, f'{db}.info.gz', check_empty=True)
		taxonlist = taxid[arg.tax] if arg.tax else None
		cls._align(query, db, outfile, taxonlist=taxid[arg.tax])

		aln_filter = cls._filter(outfile, query_cover=qc, subject_cover=sc, evalue=evalue, top=1)

		cls._anno(db+'.info.gz', aln_filter)

	def index_diamond_db(self, db, taxid:bool=True):
		"""Index the reference with diamond.
		:params
		db : str
			Input fasta to index as database.
		"""

		cmd = ['diamond', 'makedb', '--in', db, '-d', db]
		taxs = cls._get_tax_files()
		if taxs:
			cmd += ['--taxonmap', taxs[0], '--taxonnodes', taxs[1], '--taxonnames', taxs[1]]
		logger.info("*INDEXING* a db is *RUNNING*. Command is %s.", " ".join(cmd))
		path.exe_proc(cmd, shell=False)
		logger.info("*INDEXING* is *FINISHED*")

	def _align(self, query, db, aln_out, cpus=10, taxonlist=None):
		"""
		Performing alignment with diamond.

		Parameters
		----------
		query : str
			query for alignment
		db : str
			db for alignment
		aln_out : str
			output file of alignment
		taxonlist : int(taxid)
			taxonlist that align to, can be list. default None.

		Result
		------
			Output diamond alignment tabular format result.
		"""

		aln_out = path.real_path(aln_out)
		files.check_exist(db + '.info.gz', check_empty=True)
		out_columns = ['qseqid', 'qlen', 'qcovhsp',
						'sseqid', 'slen', 'scovhsp',
						'pident', 'length', 'mismatch',
						'gapopen', 'qstart', 'qend',
						'sstart', 'send', 'evalue', 'bitscore']

## no --compress, do not compress the align result this time. Jie.
		cmd = ['diamond', 'blastp', '--query', query, '--db', db,
				'--sensitive',
				'-k', '5', '-e', '0.00001', '--id', '30',
				'--tmpdir', path.get_path(aln_out), '--log',
				'-o', aln_out, '-p', str(cpu)]
		cmd += ['--outfmt', '6'] + out_columns
		if taxonlist:
			cmd += ['--taxonlist', taxonlist]

		logger.info("*ALIGNING* querys is *RUNNING*. Comannd is %s", " ".join(cmd))
		path.exe_proc(cmd, shell=False)
		files.check_empty(aln_out)

		logger.info("*ALIGNING* is *FINISHED*.")


	def _filter(self, aln, query_cover=50, subject_cover=50, evalue=1e-10, top=1):
		"""
		Filter alignments according to the identity, evalue, query_cover and subject_cover.

		Parameters
		----------
		aln : str
			tabular alignments file.
		identity : float
			minimum identity to keep an alignment, default is 30
		query_cover : float
			minimum percentage of align length cover the query, default 50.
		subject_cover : float
			minimum percentage of a lign length cover the subject, default 50.
		evalue : float
			minimum evalue to keep an alignment, default 1e-10
		top : int
			number of hits of each query to keep, default is 1.

		Results
		-------
			generate a file name suffix with *.filter
		"""

		columns = ['qseqid', 'qlen', 'qcovhsp',
					'sseqid', 'slen', 'scovhsp',
					'pident', 'length', 'mismatch',
					'gapopen', 'qstart', 'qend',
					'sstart', 'send', 'evalue', 'bitscore']

		logger.info("*FILTERING* alignments is *RUNNING*.")
		aln_filter = files.get_prefix(aln, include_path=True) + '.filter'
		aln = pd.read_csv(aln, sep='\t', header=None, index_col=None)
		aln.columns = columns

		all_index = unique(aln.qseqid, return_counts=True)
		aln_top = pd.DataFrame(columns=columns)
		for idx, count in zip(all_index[0], all_index[1]):
			aln_idx = aln[aln.qseqid == idx]
			if count > top:aln_idx = aln_idx.iloc[0:top, :]
			aln_top = aln_top.append(aln_idx)

		aln_top = aln_top[aln_top.pident>=float(identity)]
		aln_top = aln_top[aln_top.evalue<=float(evalue)]

		if query_cover:
			aln_top = aln_top[aln_top.qcovhsp>=float(query_cover)]
		if subject_cover:
			aln_top = aln_top[aln_top.scovhsp>=float(subject_cover)]
		aln_top.to_csv(aln_filter, sep='\t', index=False)

		logger.info("*FILTERING* alignments is *FINISHED*.")
		return aln_filter

	def _anno(self, db_info, aln):
		"""
		Anno alignments.

		Parameters
		----------
		db_info : str
			db.info contains annotation information.
		aln_filter : str
			alingments filtered file to anno.
		"""

		logger.info("*ANNOTATING* alignments is *RUNNING*.")
		aln_file = aln + '.anno.xls'
		files.check_exist(db_info, check_empty=True)
		aln = pd.read_csv(aln, sep="\t", header=0, index_col=0)

		db_info_reader = pd.read_csv(db_info, sep='\t', header=0, index_col=0, iterator=True)
		aligned_genes = list(aln.sseqid)
		db_info = pd.DataFrame()
		while True:
			try:
				chunk = db_info_reader.get_chunk(10000000)
			# print(chunk.index)
			# print(chunk.columns)
			# print([i for i in chunk.index if i in list(aln.sseqid)])
				chunk = chunk.loc[[i for i in chunk.index if i in aligned_genes]]
				db_info = db_info.append(chunk)
			except StopIteration:
				logger.info("Finished looping the db info in.")
				break
		db_info.columns = chunk.columns

		aln_anno = aln.merge(db_info, left_on='sseqid', right_index=True)
		aln_anno.to_csv(aln_file + '.anno.xls', sep='\t')
		db_info_simple = db_info[['Protein names', 'Gene names', 'EC number', 'Cross-reference (CAZy)', 'Cross-reference (KO)']]
		aln_anno_simple = aln[['sseqid', 'pident']].merge(db_info_simple, left_on='sseqid', right_index=True)
		aln_anno_simple.columns = ['sseqid', 'pident', 'Protein names(%s)'%d_type, 'Gene names(%s)'%d_type, 'EC(%s)'%d_type, 'CAZy(%s)'%d_type, 'KO(%s)'%d_type]

		aln_anno_simple['KO(sprot)'] = [','.join(i.strip(';').split(';')) if type(i) is str else i for i in aln_anno_simple['KO(sprot)']]

		aln_anno_simple.to_csv(aln_file + '.anno.xls.simple.xls', sep="\t")
		logger.info("*ANNOTATING* is *FINISHED*.")

	def _get_tax_files(self):
		"""prot.accession2taxid.gz, nodes.dmp and names.dmp to index db."""
		db = os.environ.get("DIAMOND_MAKEDB_PATH")
		path.check_exist(db, check_empty=True)
		if db:
			taxonmap = db+'/prot.accession2taxid.gz'
			taxonnodes = db+'/taxdmp/nodes.dmp'
			taxonnames = db+'/taxdmp/names.dmp'
			path.check_file_exist([taxonmap, taxonnodes, taxonnames])
			return [taxonmap, taxonnodes, taxonnames]
		else:
			return None

if __name__ == '__main__':
	import sys
	if len(sys.argv) == 1:
		sys.exit(sys.argv[0] + ' [query] [ref] [out] [taxon] [sprot/trembl]')

	query = sys.argv[1]
	ref = sys.argv[2]
	out = os.path.realpath(sys.argv[3])
	taxon = sys.argv[4]
	d_type= sys.argv[5]
	cpus = 10

	if not os.path.isfile(ref + '.dmnd'):
		diamond.index(ref)

	diamond.align(query, ref, out, taxon)
	diamond.sift(out, query_cover=50, subject_cover=50, evalue=1e-10)
	diamond.anno(ref+'.info.gz', out+'.filter', d_type)
