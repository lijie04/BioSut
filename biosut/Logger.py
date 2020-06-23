#####################################################################################
#																					#
# log.py - set up globally log file													#
#																					#
#####################################################################################

import logging

logger = logging.getLogger(__name__)

class log(object):

	def __init__(self, logname, logger=None)
		
		self.logname = logname
		self.logger = logging.getLogger(__name__)
		self.logger.setLevel(logging.DEBUG)

		# make a handle, write log
		fh = logging.FileHandler(self.logname, 'a', encoding='utf-8')
		fh.setLevel(logging.INFO)

		ch = logging.StreamHandler()
		ch.setLevel(logging.INFO)
	
		# output format
		formatter = logging.Formatter('[%(asctime)s] %(filename)s->%(funcName)s line:%(lineno)d [%(levelname)s]%(message)s')
		fh.setFormatter(formatter)
		ch.setFormatter(formatter)

		self.logger.addHandler(fh)
		self.logger.addHandler(ch)
	
		fh.close()
		ch.close()

	def getlog(self):
		return self.logger


