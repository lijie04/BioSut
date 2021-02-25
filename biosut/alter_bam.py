"""
The :mod:`biosut.alter_bam` includes some bam operation.
"""

# Author: Jie Li <mm.jlli6t@gmail.com>
# License: GNU v3.0
# Copyrigth: 2015 -


from re import findall

from . import gt_exe
from . import alter_seq
from . import io_seq

def recover_reads(bam, ref, out_prefix, **kargs):
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

    def _parse_samtools_info(err):
#       	return findall('discarded (\d+) singletons\\n.*processed (\d+) reads', err.decode())[0]
        stat = {}
#       	stats['singletons'] = {}
#   	stats['retreived reads'] = {}
        flags = ['paired', 'forward', 'reverse']

        for flg, err in zip(flags, errs):
            err = findall('discarded (\d+) singletons\\n.*processed (\d+) reads', err.decode())[0]
            #stats['singletons'][flg] = err[0]
            #stats['retreived reads'][flg] = err[1]
            stats[flg] = [err[0], err[1]]
        return stat

    flag = {
            'secondary': ' -F 0x100 ',
            'qcfail': ' -F 0x200 ',
            'duplicates': ' -F 0x400 ',
            'supplementary': ' -F 0x800 '
            }

    gt_exe.is_executable('samtools')
    ref = ' '.join(io_seq.seq_to_dict(ref).keys())
    fq1 = f'{out_prefix}.1.fq'
    fq2 = f'{out_prefix}.2.fq'
    fq_forward = f'{out_prefix}.forward.fq'
    fq_reverse = f'{out_prefix}.reverse.fq'

    cmd = 'samtools view -b -F 0x4'
    if kargs['secondary']:cmd += flag['secondary']
    if kargs['qcfail']:cmd += flag['qcfail']
    if kargs['duplicates']: cmd += flag['duplicates']
    if kargs['supplementary']:cmd += flag['supplementary']

    samtools_info = []

    cmd_pair = cmd + f'-f 0x2 {bam} {ref}|samtools fastq -1 {fq1} -2 {fq2}'
    out, err = gt_exe.exe_cmd(cmd_pair)
    samtools_info.append(err)
    alter_seq.reorder_PE_fq(fq1, fq2, outdir=None)

    cmd_forward = cmd + f'-F 0x2 -F 0x10 {bam} {ref}|\
                        samtools fastq -o {fq_forward} -0 {fq_forward}'
    out, err = gt_exe.exe_cmd(cmd_forward)
    samtools_info.append(err)

    cmd_reverse = cmd + f'-F 0x2 -f 0x10 {bam} {ref}|\
                        samtools fastq -o {fq_reverse} -0 {fq_reverse}'
    out, err = gt_exe.exe_cmd(cmd_reverse)
    samtools_info.append(err)

    cmd = f'gzip -f {fq1} {fq2} {fq_forward} {fq_reverse}'
    gt_exe.exe_cmd(cmd)

    stat = _parse_samtools_info(samtools_info)
    with open(f'{out_prefix}.recover_reads.stat.xls', 'w') as stat_out:
        stat_out.write('\tsingletons\trecover reads\n')
        for flg in stat:
            stat_out.write(f'{flg}\t{stat[flg][0]}\t{stat[flg][1]}\n')
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

def filter_bam(infile:str=None, outfile:str=None, bam:bool=False):
    """
    Filter non-interested alignments.

    Parameters
    ----------
    infile : str
        input BAM/SAM file.
    outfile : str
        output BAM/SAM file.
    bam : bool, default is False
        Input file is bam or not.
    """

def infer_identity(aln):
	"""
    Infer identity of alignment, hard-clipped bases were not included in denominator, \
    but soft-clipped bases are included, as I consider soft-clipp as gap-open

    Parameters
    ----------
    aln : str
        alingment line from sam

    Return
    ------
        Return identity value.
    """
	matched_bases = [int(i) for i in findall(r"(\d*)M", aln.cigarstring)]
    ## use infer_query_length() as query may been trimmed while mapping (not trimming while QC)
    # actually should include trimmed bases here? 2020-08-29
	return sum(matched_bases)*100 / aln.infer_query_length()

def infer_alignlen_ratio(aln):
	"""
    Infer alignment length ratio, weighting by denominator \
    has hard-clipped bases included to infer the real read length.
    """
	return aln.query_alignment_length*100 / aln.infer_read_length() ## here include hard-clip

#def judge_insert(aln):

#def judge_edge_alignment(aln):
