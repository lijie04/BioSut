"""
The :mod:`biosut.bamer` includes some bam operation,
make sure samtools is on system path.
"""

# Author: Jie Li <jlli6t at gmail.com>
# License: GPLv3.0
# Copyright: 2019 -

import os
import pysam as ps
from re import findall
from . import biosys as bs
from . import bioseq

# TODO: the whole module, try to use pysam instead, for more flexibility.
def recover_reads(bam, ref, out_prefix, output_bam: bool = True,
                  output_fq: bool = False, **kargs):
    """Recover reads from bam that mapped to specific reference.
    Args:
        bam: str
            bam file to process.
        ref: str
            reference fasta file.
        out_prefix : str
            output prefix to use for output files, in absolute path.
        output_bam: = bool, `True`.
            set to output bam files.
        output_fq: = bool, `False`.
            set to output fastq files.
        **kargs:
            secondary = bool, default `True`.
                set to keep secondary alignments.
            qcfail = bool, default `True`.
                set to keep qcfail alignments according to samtools principle.
            duplicates = bool, default `True`.
                set to keep duplicates alignments.
            supplementary = bool, default `True`.
                set to keep supplementary alignments.
            reorder = bool, default `True`.
                set to reorder paired reads, as they are not in original order.

    Returns:
        Out put extracted reads into paired 1&2 and unpaired 3 types.
        *.1.fastq.gz, *.2.fastq.gz, *.forward.fastq.gz, *.reverse.fastq.gz
        *stats.xls file to tell you how many reads you have extracted.
        if output_fq set to True, then stat dict will be returned.
    """

    def _parse_samtools_info(errs):
        stats = {}
        flags = ['paired', 'forward', 'reverse']
        for fl, e in zip(flags, errs):
            pattern = f'discarded (\d+) singletons\\n.*processed (\d+) reads'
            e = findall(pattern, e.decode())[0]
            # stats['singletons'][flag] = err[0]
            # stats['retrieved reads'][flag] = err[1]
            stats[fl] = [e[0], e[1]]
        return stat

    flag = {
        'secondary': ' -F 0x100 ',
        'qcfail': ' -F 0x200 ',
        'duplicates': ' -F 0x400 ',
        'supplementary': ' -F 0x800 '
    }

    bs.is_executable('samtools')
    ref = ' '.join(bioseq.seq2dict(ref).keys())
    fq1 = f'{out_prefix}.1.fastq'
    fq2 = f'{out_prefix}.2.fastq'
    fq_fwd = f'{out_prefix}.fwd.fastq'
    fq_rev = f'{out_prefix}.rev.fastq'
    bam_pair = f'{out_prefix}.pair.bam'
    bam_fwd = f'{out_prefix}.fwd.bam'
    bam_rev = f'{out_prefix}.rev.bam'

    cmd = 'samtools view -b -F 0x4'
    if kargs['secondary']: cmd += flag['secondary']
    if kargs['qcfail']: cmd += flag['qcfail']
    if kargs['duplicates']: cmd += flag['duplicates']
    if kargs['supplementary']: cmd += flag['supplementary']

    samtools_info = []
    # cmd_pair = cmd + f'-f 0x2 {bam} {ref}|samtools fastq -1 {fq1} -2 {fq2}'
    # split cmd to generate bam and fastq files simultaneously, Jie
    cmd_pair = cmd + f'-f 0x2 {bam} {ref} -o {bam_pair}'
    bs.exe_cmd(cmd_pair, shell=True)
    if output_fq:
        cmd_pair = f'samtools fastq -1 {fq1} -2 {fq2} {bam_pair}'
        # err info is from samtools fastq pipe
        out, err = bs.exe_cmd(cmd_pair, shell=True)
        samtools_info.append(err)

        # because recovered reads in different order?
        # TODO: need to manually check the order, if this command is necessary.
        bioseq.sort_pe_fq(fq1, fq2, outdir=None)

    # cmd_fwd = cmd + f'-F 0x2 -F 0x10 {bam} {ref}|' \
    #                f'samtools fastq -o {fq_fwd} -0 {fq_fwd}'
    cmd_fwd = cmd + f'-F 0x2 -F 0x10 {bam} {ref} -o {bam_fwd}'
    bs.exe_cmd(cmd_fwd, shell=True)
    if output_fq:
        cmd_fwd = f'samtools fastq -o {fq_fwd} -0 {fq_fwd} {bam_fwd}'
        out, err = bs.exe_cmd(cmd_fwd, shell=True)
        samtools_info.append(err)

    # cmd_rev = cmd + f'-F 0x2 -f 0x10 {bam} {ref}|' \
    #                f'samtools fastq -o {fq_reverse} -0 {fq_reverse}'
    cmd_rev = cmd + f'-F 0x2 -f 0x10 {bam} {ref} -o {bam_rev}'
    bs.exe_cmd(cmd_rev, shell=True)
    if output_fq:
        cmd_rev = f'samtools fastq -o {fq_rev} -0 {fq_rev} {bam_rev}'
        out, err = bs.exe_cmd(cmd_rev, shell=True)
        samtools_info.append(err)

    if not output_bam:
        cmd = f'rm {bam_pair} {bam_fwd} {bam_rev}'
        bs.exe_cmd(cmd, shell=True)

    if output_fq:
        cmd = f'gzip -f {fq1} {fq2} {fq_fwd} {fq_rev}'
        bs.exe_cmd(cmd)

        stat = _parse_samtools_info(samtools_info)
        with open(f'{out_prefix}.recover_reads.stat.xls', 'w') as stat_out:
            stat_out.write('\tsingletons\trecover reads\n')
            for flg in stat:
                stat_out.write(f'{flg}\t{stat[flg][0]}\t{stat[flg][1]}\n')
        return stat


def sam2bam(sam, sort: bool = True, overlay: bool = False):
    """
    Convert sam to bam file and/or sort it.
    Args:
        sam: FILE
            input sam file.
        sort: bool, default `True`
            set to sort the bam file.
        overlay: bool, default `False`
            set to overlay existed result, if there is any.
    Returns:
        Bam or sorted bam file.
    """
    prefix = bs.remove_suffix(sam, include_path=True)
    cmd_srt = f'samtools view -@ 10 -bS {sam} |' \
              f'samtools sort -@ 10 - > {prefix}.srt.bam'
    cmd_bam = f'samtools view -@ 10 -bS {sam} -o {prefix}.bam'

    if overlay:
        bs.exe_cmd(cmd_srt, shell=True)
        return f'{prefix}.srt.bam'

    if sort:
        if not os.path.isfile(f'{prefix}.srt.bam'):
            bs.exe_cmd(cmd_srt, shell=True)
        return f'{prefix}.srt.bam'
    else:
        if not os.path.isfile(f'{prefix}.bam'):
            bs.exe_cmd(cmd_bam, shell=True)
        return f'{prefix}.bam'


def is_sorted(bam):
    """
    Check whether a bam is sorted.
    Args:
        bam: FILE
            input bam file.

    Returns:
        Bool value
    """
    cmd = f'samtools view -H {bam} | head -n 1'
    out, err = bs.exe_cmd(cmd)
    if 'unsorted' in out.decode(): return False
    return True


# code relies on samtools, so please add samtools in.
def sort_bam(bam, overlay: bool = False):
    """
    Sort input bam file.
    Args:
        bam: FILE
            input bam file.
        overlay: bool, default `False`.
            set to overlay existed result.

    Returns:
        sorted bam file.
    """
    prefix = bs.remove_suffix(bam)
    cmd = f'samtools sort -@ 10 {bam} -o {prefix}.srt.bam'
    if overlay or not os.path.isfile(f'{prefix}.srt.bam'):
        bs.exe_cmd(cmd, shell=True)
        # return f'{prefix}.srt.bam'
    return f'{prefix}.srt.bam'


def index_bam(bam, overlay: bool = False):
    """
    Index bam file. Will sort the bam first if it is unsorted.
    Args:
        bam: = FILE
            input bam file.
        overlay: boolean, default `False`
            set to overlay existed bam file.
    Returns:
        indexed bam file.
    """
    def _index_bam(bam_):
        cmd = f'samtools index -b -@ 10 {bam_}'
        bs.exe_cmd(cmd, shell=True)

    if overlay or not os.path.isfile(f'{bam}.bai'):
        if is_sorted(bam):
            _index_bam(bam)
            return bam
        srt_bam = sort_bam(bam, overlay=True)
        _index_bam(srt_bam)
        return srt_bam

    if os.path.isfile(f'{bam}.bai'): return bam

#    try:
#        _index_bam(bam)
#        return bam
#    except FileNotFoundError:  # samtools raises this error when bam is unsort.
#        srt_bam = sort_bam(bam)
#        _index_bam(srt_bam)
#        return srt_bam


def index_bam_using_pysam(bam, overlay: bool = False):
    """
    Index bam file. If the input bam file is not sorted, will sort it first.
    Args:
        bam: FILE
            input bam file.
        overlay: boolean, default `False`
            set to overlay existed bam file.

    Returns:
        indexed bam file.
    """

    if overlay:
        if is_sorted(bam):
            ps.index(bam)
            return bam
        srt_bam = sort_bam(bam, overlay=True)
        ps.index(srt_bam)
        return srt_bam

    if os.path.isfile(f'{bam}.bai'): return bam

    try:
        ps.index(bam)
        return bam
    except ps.utils.SamtoolsError:
        srt_bam = sort_bam(bam)
        ps.index(srt_bam)
        return srt_bam
