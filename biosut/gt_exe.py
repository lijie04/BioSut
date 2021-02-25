"""
The :mod:`biosut.gt_exe` includes functions relate to command execution.
"""

# Author: Jie Li (mm.jlli6t@gmail.com)
# License: GNU v3.0
# Copyrigth: 2015 -

import os
import sys
import subprocess as sp
from loguru import logger

def is_executable(*prog):
	"""
	Check whether program (s) exists in system path.

	Paramters
	---------
	prog : str
		Program (s) that will be check.

	Result
	------
		Exit if program (s) is not exist in system path.
	"""
	for p in prog:
		code = sp.run(['which', p], stdout=sp.PIPE, stderr=sp.STDOUT).returncode

		if code:
			logger.error('Program * {p} * is not found')
			sys.exit()

def exe_cmd(cmd, shell:bool=True):
	"""
	Executing your command.

	Parameters
	----------
	cmd : str, or list
		Command will be executed.
	shell : bool, default True
		Set to False while cmd is a list.

	Results
	-------
		Execute command and return output result and error messages.
	"""

	proc = sp.Popen(cmd, shell=shell, stdout=sp.PIPE, stderr=sp.PIPE)
	out, err = proc.communicate()
	if proc.returncode != 0:
		logger.error('Error encountered while executing:{cmd}\nError message:{err}')
		sys.exit()
	return out, err
