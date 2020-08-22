"""
The :mod:`biosut.logger` for logger related functions.
"""

# Author: Jie Li (mm.jlli6t@gmail.com)
# License: GNU v3.0

import logging

class creat_logger(object):

    def __init__(self, logger=None, logfile=None):
        """
        Setup logger file.
        """
        self.logger = logging.getLogger(logger)
        self.logfile = logfile
        self.logger.setLevel(logging.DEBUG)

        log_fh = logging.FileHandler(self.logfile, 'a', encoding='utf-8')
        log_fh.setLevel(logging.INFO)

        log_cs = logging.StreamHandler()
        log_cs.setLevel(logging.INFO)

        formatter = logging.Formatter('[%(asctime)s] %(filename)s->%(funcName)s->line:%(lineno)d [%(levelnames)]%(message)s')
        log_fh.setFormatter(formatter)
        log_cs.setFormatter(formatter)
