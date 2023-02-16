"""
The :mod:`biosut.bioseq` includes utilities to operate sequence files.
"""

# Author: Jie Li <jlli6t at gmail.com>
# License: GPLv3.0
# Copyright: 2018-

import os
import sys
import gzip
from loguru import logger
import pandas as pd

from re import findall, match
from .biosys import open_file, remove_suffix, check_file

# copy-and-paste from https://github.com/lh3/readfq/blob/master/readfq.py
# add chop_comment part by jlli6t
def iterator(fh, chop_comment: bool = False):
    """
    Sequence iterator.
    Args:
        fh: FILE HANDLE
            input file handle of sequence file.
        chop_comment: boolean, default `False`
            set to chop comment in sequence id.

    Returns:
        Generates a iterator.
    """
    last = None  # this is a buffer keeping the last unprocessed line
    while True:  # mimic closure; is it a bad idea?
        if not last:  # the first record or a record following a fastq
            for l in fh:  # search for the start of the next record
                if l[0] in '>@':  # fasta/q header line
                    last = l[:-1]  # save this line
                    break
        if not last: break
        # name, seqs, last = last[1:], [], None  
        # keep comment of seq id. Jie
        name = last[1:].partition(' ')[0] if chop_comment else last[1:]
        seqs, last = [], None
        for l in fh:  # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+':  # this is a fasta record
            yield name, ''.join(seqs), None  # yield a fasta record
            if not last: break
        else:  # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fh:  # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs)  # yield a fastq record
                    break
            if last:  # reach EOF before reading enough quality
                yield name, seq, None  # yield a fasta record instead
                break

def fq2fa(infq, out_fa: str = 'test'):
    """
    Convert FASTQ to FASTA format.
    Args:
        infq: FILE
            input FASTQ (.gz) file.
        out_fa: FILE
            output FASTA file.

    Results:
        Generate a FASTA file.
    """
    check_file(infq, check_empty=True)
    fh = open_file(infq)
    with open(out_fa, 'w') as outf:
        for t, seq, _ in iterator(fh):
            outf.write(f'>{t}\n{seq}\n')
    fh.close()

# TODO: test the code
def count_symbol(string, symbol:str = 'GC', ignore_case: bool = True, 
                seperate_symbol: bool = True):
    """
    Count specified symbol in a string.
    Args:
        in_string: str
            input string. Usually it is a sequence.
        symbol: str, default is `GC`
            symbol to count in string.
        ignore_case: boolean, default is `True`
            set to ignore case of letters, both in string and symbol.
        seperate_symbol: boolean, default is `True`
            set to seperate specified symbol when count. 
            e.g. When specified symbol is 'GC', 
            the func will return number of G+C if `seperate_symbol=True`, 
            otherwise return number of GC.

    Returns:
        number of specified symbol in a given string.
    """
    if ignore_case: string, symbol = string.upper(), symbol.upper()
    count = 0
    if seperate_symbol:
        symbol = list(symbol)
        if len(set(symbol)) < len(symbol):
            logger.error(f'Repeat symbols in {symbol}.')
            sys.exit()
        count = 0
        for s in symbol:
            count += string.count(s)
        return count
    else:
        return string.count(symbol)

# TODO: test the code and also the SPEED.
def seq2dict(inseq: str, sequence: bool = True, length: bool = True, 
            gc: bool = True, qual: bool = False, dataframe: bool = False):
    """
    Read sequences file into dict with sequence id as key and specified value.
    Args:
        inseq: FILE
            input FASTA/FASTQ (.gz) file.
        sequence: bool, default is `True`
            set to output sequence.
        length: boolean, default is 'True'
            set to output length of each sequence.
        gc: boolean, default is `True`
            set to output gc of each sequence.
        qual: boolean, default is `False`.
            set to output qual,
        dataframe: boolean, default is `False`.
            set to convert dict to dataframe before return.

    Returns:
        a dataframe with sequence id as the first column.
    """
    seq_dict = {}
    check_file(inseq, check_empty=True)
    fh = open_file(inseq)
    if sequence:seq_dict['sequence'] = {}
    if length:seq_dict['length'] = {}
    if gc:seq_dict['gc'] = {}
    if qual:seq_dict['qual'] = {}

    for t, seq, _ in iterator(fh):
        if sequence:seq_dict['sequence'][t] = seq
        if length:seq_dict['length'][t] = len(seq)
        if gc:seq_dict['gc'][t] = count_symbol(seq, symbol='GC', seperate_symbol=True)
        if qual:seq_dict['qual'][t] = _
    fh.close()
    if dataframe:return pd.DataFrame.from_dict(seq_dict)
    return seq_dict

def assess_genome(genome, min_len: int = 500):
    """
    Assess genome and return genome traits.
    Args:
        genome: FILE.
            input a genome FASTA file.
        min_len:= int, default `500`.
            minimum length to retain a sequence.

    Returns:
        return a list contains [genome size, contig number, n50,
        maximal contig length, minimal contig length, gap number,
        N number, gc ratio].
    """
    check_file(genome)
    fh = open_file(genome)
    n, gap, gc, ctg_num, ctg_len, n50 = 0, 0, 0, 0, [], None
    for t, seq, _ in iterator(fh):
        if len(seq) < min_len: continue
        ctg_num += 1
        ctg_len.append(len(seq))
        gap += len(findall('N+', seq))  # to test the findall function.
        n += seq.count('N')
        gc += count_symbol(seq, symbol='GC', seperate_symbol=True)

    genome_size = sum(ctg_len)
    gc = round(gc/genome_size*100., 2)

    ctg_len.sort(reverse=True)
    sum_len = 0
    for leng in ctg_len:
        sum_len += leng
        while sum_len >= genome_size*0.5: n50 = leng
    return [genome_size, ctg_num, n50, max(ctg_len), min(ctg_len), gap, n, gc]


# TODO: to check the code is necessory
def select_seq_len(inseq=None, outseq=None,
                   max_len: int = 0, min_len: int = 0,
                   rm_longest_percent: float = 0,
                   rm_shortest_percent: float = 0,
                   output_qual: bool = False):
    """
    Select sequences according the length.
    Args:
        inseq: FILE
            input FASTA/Q (.gz) file
        outseq: FILE
            output sequence file name.
        max_len:= int, default `0`.
            maximal length of sequence to keep, 0 means keep all.
        min_len:= int, default `0`.
            minimal length of sequence to keep, 0 means keep all.
        rm_longest_percent:= float, default `0.0`.
            longest n% sequences will be excluded. 0 means keep all.
        rm_shortest_percent:= float, default `0.0`.
            shortest n% sequences will be excluded. 0 means keep all.
        output_qual:= bool, default `False`.
            set to output quality line as well.

    Results:
        Generates number trimmed FASTA/Q file.
    """
    def cal_percent(inseq_, min_len_, longest_perc, shortest_perc):
        all_seq_len = list(seq2dict(inseq_, min_len=min_len_).values())
        all_seq_len.sort(reverse=True)
        seq_num = len(all_seq_len)
        long_ = all_seq_len[round(longest_perc * seq_num + 0.5)]
        short_ = all_seq_len[seq_num - round(shortest_perc * seq_num + 0.5)-1]
        return long_, short_

    if not outseq:
        logger.error(f'U have to set outseq file.')
        sys.exit()

    if rm_longest_percent or rm_shortest_percent:
        long, short = cal_percent(inseq, min_len, rm_longest_percent,
                                  rm_shortest_percent)

    fh = open_file(inseq)
    with open(outseq, 'w') as outf:
        for t, seq, qual in iterator(fh):
            if max_len and len(seq) > max_len: continue
            if min_len and len(seq) < min_len: continue
            if rm_longest_percent and len(seq) > long: continue
            if rm_shortest_percent and len(seq) < short: continue
            if output_qual:
                outf.write(f'@{t}\n{seq}\n+\n{qual}\n')
            else:
                outf.write(f'>{t}\n{seq}\n')
    fh.close()

# TODO: test the code
def split_seq(in_fa: str, out_fa: str, symbol: str = 'N', exact: bool = True):
    """
    Split sequences according to specified symbol. Such as Ns.
    Args:
        in_fa: FILE
            input FASTA file.
        out_fa: FILE
            output FASTA file name.
        symbol:= str, default `N`.
            symbol for using to fragment sequence.
        exact:= bool, default `True`.
            set to indicting exact symbol, otherwise, will be fuzzy matching.
            e.g, set symbol to NN, and exact=True,
            program will not recognize NNN or NNNN as a split site.

    Results:
        Output split sequences into specified out_fa FILE.
    """
    symbol_len = len(symbol)
    symbol += '+'  # make a 're' match to indicate one or more symbol
    fh = open_file(in_fa)
    # out = open(out_fa, 'w')
    print(symbol, symbol_len)
    with open(out_fa, 'w') as out:
        for t, seq, _ in iterator(fh):
            c, start, end = 0, 0, 0
            gaps = findall(symbol, seq)

            if len(gaps) == 0:
                out.write(f'>{t}\n{seq}\n')
                continue

            for gap in gaps:
                pos = seq[end:].find(gap)
                end += pos
                # use symbol_len to replace, judge whether to stop here or not.
                if len(gap) == symbol_len:  # exact a 'gap' to split fasta
                    c += 1
                    out.write(f'>{t}_{c}|len={end-start}\n{seq[start:end]}\n')
                    start = end + len(gap)
                    end += len(gap)
                    continue
                if exact:  # number os `symbol` is more than expected.
                    end += len(gap)
                    continue
                c += 1
                out.write(f'>{t}_{c}|len={end-start}\n{seq[start:end]}\n')
                start = end + len(gap)
                end += len(gap)
            # output the last one, as n gaps cut sequences into n+1 sequences.
            out.write(f'>{t}_{c + 1}|len={len(seq) - start}\n{seq[start:]}\n')
    fh.close()


# TODO: to check if the code is necessary.
def sort_pe_fq(infq1: str, infq2: str, outdir=None):
    """
    Order pair-end FASTQ files to make fq1 & fq2 in same order.
    Args:
        infq1: FILE
            input FASTQ 1 file (.gz)
        infq2: FILE
            input FASTQ 2 file (.gz)
        outdir:= str, default `none`
            Outdir to output files.
            With outdir not specified, old files will be overlaid.

    Results:
        Output pair-end FASTQ files that in same order.
    """
    if infq1[-3:] == '.gz': infq1 = remove_suffix(infq1)
    if infq2[-3:] == '.gz': infq2 = remove_suffix(infq2)

    if outdir:
        fq1_out = open(f'{outdir}/{os.path.basename(infq1)}', 'w')
        fq2_out = open(f'{outdir}/{os.path.basename(infq2)}', 'w')
    else:  # output fastq files will replace the old ones.
        fq1_out = open(infq1, 'w')
        fq2_out = open(infq2, 'w')

    # TODO: fq1 and fq2 should be in same length?
    #  if not, then it is the problem of samtools? Jie, 2021-10-30
    fq1 = seq2dict(infq1, qual=True)
    fq2 = seq2dict(infq2, qual=True)

    # print('Reordering fastq sequences id.')
    for sid in fq1.keys():
        fq1_seq_qual = fq1[sid].split(',')
        fq2_seq_qual = fq2[sid].split(',')
        fq1_out.write(f'@\n{sid}\n{fq1_seq_qual[0]}\n+\n{fq1_seq_qual[1]}\n')
        fq2_out.write(f'@\n{sid}\n{fq2_seq_qual[0]}\n+\n{fq2_seq_qual[1]}\n')
    fq1_out.close()
    fq2_out.close()


def extract_seq(inseq: str = None, idlist: str = None, match_out: str = None,
                outqual: bool = False, unmatch_out: str = None,
                exact_match: bool = True):
    """
    Extract sequences using sequence id specified in idlist.
    Args:
        inseq: FILE
            input FASTA/Q (.gz) file.
        idlist:= list
            idlist for sequence extracting.
        match_out:= FILE
            file name to output matched sequences.
        outqual:= boolean, default `False`.
            set to include quality line in output if input is in FASTQ format.
        unmatch_out:= FILE
            if specified, negative match sequences will be wrote into this file.
        exact_match:= boolean, default `True`
            set to do exact matching of id. Otherwise do fuzzy matching.

    Results:
        Output matched and/or negative matched sequence file.
    """
    if not match_out and not unmatch_out:
        logger.error(f'No match_out or unmatch_out file specified. Exiting...')
        sys.exit()

    if match_out: match_out_file = open(match_out, 'w')
    if unmatch_out: unmatch_out_file = open(unmatch_out, 'w')

    fh = open_file(inseq)
    for t, seq, _ in iterator(fh):
        line = outqual and f'@{t}\n{seq}\n+\n{_}\n' or f'>{t}\n{seq}\n'
        if exact_match:
            if t in idlist:
                if match_out: match_out_file.write(line)
            else:
                if unmatch_out: unmatch_out_file.write(line)
            continue

        flag = 0
        for i in idlist:
            if match(i, t):
                flag = 1
                break

        if flag and match_out: match_out_file.write(line)
        if flag and unmatch_out: unmatch_out_file.write(line)
    fh.close()
    if match_out: match_out_file.close()
    if unmatch_out: unmatch_out_file.close()


def trim_headn(inseq: str = None, outseq: str = None, outqual: bool = False):
    """
    Trim N from head of sequence.(Because sequence has N from the start?)
    Args:
        inseq:= FILE
            input FASTA/Q (.gz) file.
        outseq: FILE
            file name to output sequences,
            output will be compressed if ends with .gz
        outqual:= boolean, default `False`.
            set to output quality line if input is in FASTQ format.

    Results:
        output N-trimmed sequences into output file.
    """
    def remove_first_n(string):
        while string[0] in 'Nn':
            string = string[1:]
        return string

    if not inseq:
        logger.info('Must specify an input seq file for trim_headn.')
        sys.exit()

    if not outseq:
        logger.info('Must specify an output file for trim_headn.')
        sys.exit()

    fh = open_file(inseq)
    ofh = '.gz' in outseq and gzip.open(outseq, 'wb') or open(outseq, 'w')
    for t, seq, _ in iterator(fh):
        seq = remove_first_n(seq)
        line = f'>{t}\n{seq}\n'
        if outqual: line = f'@{t}\n{seq}\n+\n{_[len(_)-len(seq):]}\n'
        ofh.write('.gz' in outseq and line.encode() or line)
    fh.close()
    ofh.close()
