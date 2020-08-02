"""
The :mod:`biosut.alter_bam` includes some bam operation.
"""

# Author: Jie Li <mm.jlli6t@gmail.com>
# License: GNU v3.0

from re import findall

import .go_exe
import .alter_seq

def extract_reads(bam, ref, out_prefix, **kargs):
    """
    Extract reads from bam that mapped reference.

    Parameters
    ----------
    bam : str
        bam file to process
    ref : str
        referece fasta file to use
    out_prefix : str
        output prefix to use for output files, must include paths.

    secondary : bool, default is True
        keep secondary alignments.
    qcfail : bool, default is True
        set to keep qcfail alignments according to samtools principle.
    duplicates : bool, deafult is True
        set to keep duplicates alignments.
    supplementary : bool, default is True
        set to keep supplementary alignments.
    reorder : bool, default is True
        set to reorder paired sequences, as they are not in same order originally.

    Results
    -------
        Out put extracted reads into paired 1&2 and unpaired 3 types.
        *.1.fastq.gz, *.2.fastq.gz, *.forward.fastq.gz, *.reverse.fastq.gz\
        *stats.xls file to tell you how many reads you have extracted.
        stat dict will be returned, fq files won't
        ### And return paired fastq 1 &2, unpaired forward fastq 1 and reverse fastq 2 file.
    """

    flag = {
            'secondary': '-F 0X100 ',
            'qcfail': '-F 0X200 ',
            'duplicates': '-F 0X400 ',
            'supplementary': '-F 0X800 '
            }

    ref = ' '.join(io_seq.seq_to_dict(ref).keys())
    fq1 = '%s.1.fastq' % out_prefix
    fq2 = '%s.2.fastq' % out_prefix
    fq_forward = '%s.forward.fastq' % out_prefix
    fq_reverse = '%s.reverse.fastq' % out_prefix

    cmd = 'samtools view -b '
    if kargs['secondary']:cmd += flag['secondary']
    if kargs['qcfail']:cmd += flag['qcfail']
    if kargs['duplicates']: cmd += flag['duplicates']
    if kargs['supplementary']:cmd += flag['supplementary']

    samtools_info = []

    cmd_pair = cmd + '-f 0X2 %s %s|samtools fastq -1 %s -2 %s' % \
                (bam, ref, fq1, fq2)
    out, err = go_exe.exe_cmd(cmd_pair)
    samtools_info.append(err)
    alter_seq.reorder_PE_fq(fq1, fq2, outdir=None)

    cmd_forward = cmd + '-F 0X2 -F OX10 %s %s|samtools fastq -o %s -0 %s' % \
                    (bam, ref, fq_forward, fq_forward)
    out, err = go_exe.exe_cmd(cmd_forward)
    samtools_info.append(err)

    cmd_reverse = cmd + '-F 0X2 -f 0X10 %s %s|samtools fastq -o %s -0 %s' % \
                    (bam, ref, fq_reverse, fq_reverse)
    out, err = go_exe.exe_cmd(cmd_reverse)
    samtools_info.append(err)

    cmd = 'gzip -f %s -2 %s %s %s' % (fq1, fq2, fq_forward, fq_reverse)
    go_exe.exe_cmd(cmd)

    stat = _parse_samtools_info(samtools_info)
    with open(out_prefix+'.retreived_reads.stat.xls', 'w') as stat_out:
        stat_out.write('\tsingletons\tretreived reads\n')
        for flg in stat:
            stat_out.write('%s\t%s\t%s\n' % (flg, stat[flg][0], stat[flg][1]))
    return stat

def _parse_samtools_info(err):
#	return findall('discarded (\d+) singletons\\n.*processed (\d+) reads', err.decode())[0]
    stat = {}
#	stats['singletons'] = {}
#	stats['retreived reads'] = {}
    flags = ['paired', 'forward', 'reverse']

    for flg, err in zip(flags, errs):
        err = findall('discarded (\d+) singletons\\n.*processed (\d+) reads', err.decode())[0]
        #stats['singletons'][flg] = err[0]
        #stats['retreived reads'][flg] = err[1]
        stats[flg] = [err[0], err[1]]
    return stat

# deprecated, use alter_seq.reorder_PE_fq(infq1, infq2, outdir=None) instead.
#	def _re_order_pair_fq(fq1, fq2):
#		fq1_out = fq1
#		fq2_out = fq2
#		fq1 = sequtil.read_seq_to_dict(fq1, qual=True)
#		fq2 = sequtil.read_seq_to_dict(fq2, qual=True)
#	print('Reordering fastq sequences id.')
#		with open(fq1_out, 'w') as fq1_out, open(fq2_out, 'w') as fq2_out:
#			for sid in fq1.keys():
#				fq1_out.write('@\n%s\n%s\n+\n%s\n'%(sid, fq1[sid][0], fq1[sid][1])
#				fq2_out.write('@\n%s\n%s\n+\n%s\n'%(sid, fq2[sid][0], fq2[sid][1])


#def infer_identity(aln):
#	""" Infer identity of alignment, denominator not include hard-clipped bases, soft-clipped bases are included, I consider soft-clipp as gap-open"""
#	matched_bases = [int(i) for i in findall(r"(\d*)M", aln.cigarstring)] ## M indicates matched, so extracted all "Matched" bases
#	return sum(matched_bases)*100 / aln.infer_query_length()  ## use infer_query_length() as query may been trimmed while mapping (not trimming while QC)


#def infer_alignlen_ratio(aln):
#	"""Infer alignment length ratio weighting by Denominator has hard-clipped bases included to infer the real read length"""
#	return aln.query_alignment_length*100 / aln.infer_read_length() ## here I have hard-clipped bases included

#def judge_insert(aln):

#def judge_edge_alignment(aln):
