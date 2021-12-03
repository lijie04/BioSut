"""
The :mod:`biosut.mapper` includes utilities to operate sequence files.
"""

# Author: Jie Li <jlli6t at gmail.com>
# License: GPLv3.0
# Copyright: 2021


import os
from loguru import logger
from . import biosys as bs
from . import bamer


class Mapper:
    def __init__(self, fq1, fq2, ref, outdir, unfq=None,
                 min_insert_size: int = 0, max_insert_size: int = 1000,
                 cpu: int = 20, fmt: str = "sorted"):
        """
        Start the mapper tool set.
        Args:
            fq1: FILE
                input fq1 file.
            fq2: FILE
                input fq2 file.
            ref: FILE
                input reference file.
            outdir: str
                output directory.
            unfq: FILE, default `None`
                unpaired fq file.
            min_insert_size: = int, default `0`
                minimum insert size of library fragment.
            max_insert_size: = int, default `1000`
                maximum insert size of library fragment.
            cpu: = int, default `20`
                number of cpu to use.
            fmt: str, default `sorted`
                format of output file, [sam/bam/sorted].
        """
        self.fq1 = fq1
        self.fq2 = fq2
        self.ref = ref
        self.outdir = bs.sure_path_exist(outdir)
        self.unfq = unfq
        self.min_insert_size = min_insert_size
        self.max_insert_size = max_insert_size
        self.cpu = cpu
        self.fmt = fmt
        bs.check_file_exist(self.fq1, self.fq2, self.ref, check_empty=True)
        self.fq_prefix = bs.remove_suffix(self.fq1, seq=True)
        self.ref_prefix = bs.remove_suffix(self.ref)

    def bowtie2(self):
        """
        Run bowtie2 for pair-end reads and/or single reads.
        Returns:
            return the mappings in sorted bam format.
        """
        logger.info('Running bowtie2 mapping.')
        sub_outdir = f'{self.outdir}/bowtie2'
        bs.sure_path_exist(sub_outdir)
        if not os.path.isfile(f'{self.ref}.1.bt2'):
            logger.info(f'bowtie2 index of {self.ref} is not exist.')
            cmd = f'bowtie2-build {self.ref} {self.ref}'
            bs.exe_cmd(cmd)
        prefix = f'{sub_outdir}/{self.fq_prefix}.map.{self.ref_prefix}'
        cmd = f'bowtie2 -S {prefix}.sam ' \
              f'--phred33 --very-sensitive-local -p {self.cpu} ' \
              f'-I {self.min_insert_size} -X {self.max_insert_size} ' \
              f'-x {self.ref} -1 {self.fq1} -2 {self.fq2}'
        if self.unfq:
            cmd += f' -U {self.unfq}'
        cmd += f'2>{prefix}.stat'
        bs.exe_cmd(cmd, shell=True)
        logger.info('Finished bowtie2 running.')
        if self.fmt == 'sam': return f'{prefix}.sam'
        return bamer.sam2bam(f'{prefix}.sam')

    def bbmap(self):
        """
        Run bbmap for pair-end reads.
        Returns:
            return the mappings in sorted bam format.
        """
        logger.info('Running bbmap mapping.')
        sub_outdir = f'{self.outdir}/bbmap'
        bs.sure_path_exist(sub_outdir)
        prefix = f'{sub_outdir}/{self.fq_prefix}.map.{self.ref_prefix}'
        cmd = f'bbmap.sh ref={self.ref} k=13 ' \
              f'path={bs.get_file_path(self.ref)} ' \
              f'in={self.fq1} in2={self.fq2} ambiguous=random ' \
              f'pairlen={self.max_insert_size-2*150} qin=33 ' \
              f'out={prefix}.sam -Xmx40g threads={self.cpu} nodisk ' \
              f'2>{prefix}.stat'
        bs.exe_cmd(cmd, shell=True)
        logger.info('Finished bbmap running.')
        if self.fmt == 'sam': return f'{prefix}.sam'
        return bamer.sam2bam(f'{prefix}.sam')
