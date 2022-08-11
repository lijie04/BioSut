"""
The :mod:`biosut._diamond` integrated diamond related operations.
"""

# Author: Jie Li <jlli6t at gmail.com>
# License: GPLv3.0
# Copyright: 2021 -

import os
import pandas as pd
from numpy import unique
from loguru import logger

from biosut import biosys as bs


# TODO: to fix a whole bunch of bug.
class Aligner:
    taxid = {'archaea': '2157', 'bacteria': '2'}  # define tax id.
    # diamond_columns = ['qseqid', 'qlen', 'qcovhsp', 'sseqid', 'slen',
    #                   'scovhsp','pident', 'length', 'mismatch', 'gapopen',
    #                   'qstart', 'qend', 'sstart', 'send', 'evalue',
    #                   'bitscore']
    diamond_columns = ['qlen', 'qcovhsp', 'sseqid', 'slen', 'scovhsp', 'pident',
                       'length', 'mismatch', 'gapopen', 'qstart', 'qend',
                       'sstart', 'send', 'evalue', 'bitscore']

    def __init__(self, query, subject, outdir, db_type: str = 'aa',
                 cpu: int = 10, diamond_taxdb: str = 'DIAMOND_TAXDB',
                 tax: str = 'bacteria', identity: float = 30,
                 query_cov: int = 50, subject_cov: int = 50,
                 evalue: float = 1e-10, top: int = 1):
        """
        Start to run aligner module.
        Args:
            query: FILE
                input query sequence file.
            subject: FILE
                input subject sequence file.
            outdir: str
                output directory.
            db_type: = str, default `aa`
                database type, [aa/nt].
            cpu: = int, default `10`
                number of cpu to use.
            diamond_taxdb: str, default `DIAMOND_TAXDB`
                variable name of taxdb to index diamond database.
            tax: = str, default `bacteria`
                taxonomy of query. [bacteria/archaea]
        """
        self.query = query
        self.subject = subject
        self.outdir = bs.sure_path_exist(outdir)
        self.db_type = db_type
        self.cpu = cpu
        self.diamond_taxdb = diamond_taxdb
        self.tax = tax
        self.identity = identity
        self.query_cov = query_cov
        self.subject_cov = subject_cov
        self.evalue = evalue
        self.top = top
        bs.check_file_exist(self.query, self.subject, f'{self.subject}.info.gz')

    def diamond(self):
        """
        Run diamond alignments.

        Returns:
            return the filtered diamond alignments DataFrame.
        """
        def index_diamond_db(subject, diamond_taxdb):
            cmd_makedb = f'diamond makedb --in {subject} -d {subject}'
            taxs = get_tax_file(diamond_taxdb)
            if taxs:
                cmd_makedb += f' --taxonmap {taxs[0]} --taxonnodes {taxs[1]}' \
                              f' --taxonnames {taxs[2]}'
            logger.info('Indexing diamond database.')
            bs.exe_cmd(cmd_makedb, shell=True)
            logger.info('Finished indexing diamond database.')

        def get_tax_file(diamond_taxdb):
            tax_db = os.environ.get(diamond_taxdb)
            bs.check_path_exist(tax_db, check_empty=True)
            if tax_db:
                logger.info(f'Found DIAMOND_TAXDB')
                taxon_map = f'{tax_db}/prot.accession2taxid.gz'
                taxon_nodes = f'{tax_db}/taxdmp/nodes.dmp'
                taxon_names = f'{tax_db}/taxdmp/names.dmp'
                bs.check_file_exist(taxon_map, taxon_nodes, taxon_names)
                return [taxon_map, taxon_nodes, taxon_names]
            else:
                logger.info(f'DIAMOND_TAXDB not found, return None')
                return None

        # index this db
        if not os.path.isfile(f'{self.subject}.dmnd'):
            index_diamond_db(self.subject, self.diamond_taxdb)
        sub_outdir = f'{self.outdir}/diamond'
        logger.info('Running diamond alignment.')
        align_out_file = f'{sub_outdir}/diamond.align.out'
        cmd_align = f'diamond blastp --query {self.query} ' \
                    f'--db {self.subject} --sensitive -k 5 -e 0.00001 ' \
                    f'--id 30 --tmpdir {sub_outdir} ' \
                    f'--log -o {align_out_file} -p {self.cpu} ' \
                    f'--outfmt 6 {" ".join(self.diamond_columns)} ' \
                    f'--taxonlist {self.taxid[self.tax]}'
        bs.exe_cmd(cmd_align, shell=True)
        logger.info('Finished diamond alignment.')
        align_out = pd.read_csv(align_out_file, sep='\t', header=None,
                                index_col=0)
        align_out.columns = self.diamond_columns

        align_out.to_csv(align_out_file, sep='\t')
        align_filter = self.filter_aln(align_out)
        return align_filter, f'{align_out_file}.filter'

    # def blastp(self):

    def filter_aln(self, align_out):
        """
        Filter alignment results.
        Args:
            align_out: DataFrame
                input alignments DataFrame.

        Returns:
            Filtered alignment DataFrame.
        """
        logger.info('Filter alignments.')
        align_filter_file = f'{align_out}.filter'
        align_out = pd.read_csv(align_out, sep='\t', header=0, index_col=0)
        all_uniq_qseqid = unique(align_out.index, return_counts=True)
        align_top = pd.DataFrame(columns=align_out.columns)
        for idx, cnt in zip(all_uniq_qseqid[0], all_uniq_qseqid[1]):
            align_idx = align_out.loc[idx]
            if cnt > self.top: align_idx = align_idx.iloc[0:self.top, :]
            align_top = align_top.append(align_idx)

        align_top = align_top[align_top.pident >= self.identity]
        align_top = align_top[align_top.evalue <= self.evalue]
        align_top = align_top[align_top.qcovhsp >= self.query_cov]
        align_top = align_top[align_top.scovhsp >= self.subject_cov]
        align_top.to_csv(align_filter_file, sep='\t')
        logger.info('Finished filtering alignments.')
        return align_top


if __name__ == '__main__':
    import sys
    if len(sys.argv) == 1:
        sys.exit(sys.argv[0] + ' [query] [ref] [out] [sprot/trembl]')

    querys = sys.argv[1]
    ref = sys.argv[2]
    out = os.path.realpath(sys.argv[3])

    aligner = Aligner(querys, ref, out)
    diamond_align_filter = aligner.diamond()
