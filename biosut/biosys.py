"""
The :mod:`biosut.biosys` includes system functions.
"""

# Author: Jie Li (jlli6t at gmail.com)
# License: GPLv3.0
# Copyright: 2018 -

import os
import sys
import re
import gzip
import subprocess as sp
from loguru import logger


def is_executable(*prog):
    """
    Check whether program(s) exists in system path.
    Args:
        *prog: str, prog1, pro2, pro3, ...
            program (s) that will be checked.

    Results:
        Exit if program is not exist in system path.
    """
    for p in prog:
        run_prog = sp.run(['which', p], stdout=sp.PIPE, stderr=sp.STDOUT)

        if run_prog.returncode:
            logger.error(f'Program * {p} * is system path.')
            sys.exit()

def execute_cmd(cmd, shell: bool = True):
    """
    Executing your command.
    Args:
        cmd: str, list
            input command to be executed.
        shell: boolean, default `True`
            set when command is in string format.

    Results:
        execute command and return output result and error messages.
    """
    proc = sp.Popen(cmd, shell=shell, stdout=sp.PIPE, stderr=sp.PIPE)
    out, err = proc.communicate()
    if proc.returncode != 0:
        logger.error(f'Error encountered while executing: {cmd}\n'
                     f'Error message:{err}')
        sys.exit()
    return out, err

def check_path(*paths_in, check_empty: bool = False, mkdir: bool = False):
    """
    Check if path exists, exits if path does not exist. When `mkdir=True`, the path will be created if not exists.
    Args:
        *paths_in: str (s)
            path (s) that will be checked.
        check_empty: boolean, default `True`
            set to check if the path is empty. Exit if path is empty.
        mkdir: boolean, default `False`
            set to create the path if not exists.

    Returns:
        return full path(s).

    Example:
        >>> from biosut.biosys import sure_path_exist
        >>> check_path('/a/b/c', check_empty=True, mkdir=False)
    """
    def check_path_empty(dir_in):
        """
        Check if direcotry is empty.
        Args:
            dir_in: str
                directory to check.
        Results:
            exit and report error message if input directory is empty.
        """
        if not os.listdir(dir_in):
            logger.error(f'Directory {dir_in} is empty.')
            sys.exit()

    final_paths = []
    for pth in paths_in:
        pth = os.path.abspath(pth)
        final_paths.append(pth)
        if os.path.exists(pth):
            if check_empty:
                check_path_empty(pth)
        else:
            if mkdir:
                try:
                    os.makedirs(pth)
                except OSError as e:
                    logger.error(f'Error encountered: {e}, path {pth} cannot be created.')
                    sys.exit()
            else:
                logger.error(f'path {pth} does not exists and you didnot set mkdir=True to create it.')
                sys.exit()
    return len(final_paths) == 1 and final_paths[0] or final_paths

# TODO
def real_path(*paths):
    """
    Get real path, soft links will be converted to solid destination path.
    Args:
        *paths: str (s)
            input path (s).
    Returns:
        return real path as a string or path (s) as a list.

    Examples:
        >>> from biosut.biosys import real_path
        >>> a = './../bucket'
        >>> b = '/usr/bin/python'
        >>> c = 'aaa -> ../../aaa' # this is a soft link
        >>> result = real_path(a, b, c)
        >>> print(result)
        ['/path/to/../bucket', '/usr/bin/python', '/path/to/../../aaa']
    """
    final_paths = [os.path.realpath(pth) for pth in paths]
    return len(final_paths) == 1 and final_paths[0] or final_paths

def find_db(db_v: str):
    """
    Find the database indicating by db_b.
    Args:
        db_v: str
            database env variable.

    Returns:
        Return the database if found and not empty.
        Otherwise, program will exits.

    Examples:
        >>> from biosut.biosys import find_db
        >>> find_db('KEGG')
        '/full/path/to/../../aaa'
    """
    db = os.environ.get(db_v)
    if db:
        check_path(db, check_empty=True)
        return db
    logger.error(f'{db_v} could not be found.')
    sys.exit()

def check_file(*files_in, check_empty: bool = False):
    """
    Check if file (s) exists, exit if file not exist.
    Args:
        *files_in: str
            input file (s) to check.
        check_empty: boolean, default `False`
            set to check if file is empty
    Results:
        exit if the file is not exist or emtpy when check_empty=True.
    """
    def check_file_empty(file_in):
        """
        Check whether the file is empty.
        Args:
            file_in: str
                input file (s) to check emptiness.

        Results:
            exit if file is empty.
        """
        if not os.path.getsize(file_in):
            logger.error(f'{file_in} is empty. Exiting...')
            sys.exit()

    for fl in files_in:
        if not os.path.isfile(fl):
            logger.error(f'{fl} is not exist. Exiting...')
            sys.exit()
        if check_empty: check_file_empty(fl)

def remove_suffix(file_in: str = None, split_symbol: str = '.', 
                times: int = 1, include_path: bool = True):
    """
    Get prefix of file. e.g. file is `test.txt`, then return `test`.
    Args:
        file_in: = str
            input file, relative path or absolute path.
        split_symbol: = str, default `.`
            symbol to use to split file name.
        times: = int, default `1`.
            how many times supposed to chop the string.
        include_path: = boolean, default `True`
            set to include path in output.
        deprecated parameter:
        seq: = boolean, default `False`
            set to indicate input file is gzipped FASTA
    Returns:
        return prefix, w/o absolute path.
    """
    # check_file_exist(file_in)  # still get the prefix even file is not exist
    file_in = os.path.abspath(file_in)
    # return re.sub('.\d.fq.gz|.\d.fa.*.gz|.\d.fa.*|.\d.fq', '', file_in)
    # if seq: return re.sub('.\d.fq.gz|.\d.fa.*|.\d.fq', '', file_in) # can use para times=2 to control?-Jie-20230215
    file_in = file_in.split(split_symbol)
    file_in = split_symbol.join(file_in[0:len(file_in) - times])
    return include_path and file_in or os.path.basename(file_in)

# TODO
def get_file_path(*files_in):
    """
    Get absolute path of input and return.
    Args:
        *files_in: str
            input file (s), could be gzipped.

    Returns:
        absolute path (s) of the input file (s).
    """
    final_paths = [os.path.dirname(os.path.abspath(fl)) for fl in files_in]
    return len(final_paths) == 1 and final_paths[0] or final_paths

# TODO
def open_file(file_in: str = None):
    """
    Open a file and return a file handle.
    Args:
        file_in: = str
            input a file.

    Returns:
        return a file handle.
    """
    return '.gz' in file_in and gzip.open(file_in, 'rt') or open(file_in, 'r')

# TODO
def find_file(dr: str = None, suffix='fa'):
    """
    Find files with specified suffix.
    Args:
        dr: = str
            input direcotry for finding files.
        suffix: = `fa`

    Returns:
        a list of found files with full path.
    """
    return [dr + '/' + fl for fl in os.listdir(dr) if fl.endswith(suffix)]

# TODO
def parse_json(json_in):
    """
    Parse json file. (It really depends on how deep the json file is).
    Args:
        json_in: str
            input json format file.

    Returns:
        return a string but formatted like a dataframe.
    """
    import json
    fh = open(json_in)
    # with open(json_in) as fp:
    # out.write('Head line.\n')
    json_file = json.load(fh)
    # out_df = pd.DataFrame(index=None, columns=None)
    out_st_df = ''
    n0 = json_file['name']
    for cld1 in json_file['children']:
        n1 = cld1['name']
        for cld2 in cld1['children']:
            n2 = cld2['name']
            for cld3 in cld2['children']:
                n3 = cld3['name']  # "ko01000" is strange. what this for?
                if 'children' not in cld3:
                    # out.write(f'{n0}\t{n1}\t{n2}\t{n3}\t{n4}\n')
                    # out.write(f'{n0}\t{n1}\t{n2}\t{n3}\n')
                    # out_df.append(pd.DataFrame([n0, n1, n2, n3]).T)
                    out_st_df += f'{n0}\t{n1}\t{n2}\t{n3}\n'
                else:
                    # print(cld3)
                    for cld4 in cld3['children']:
                        n4 = cld4['name']
                        if 'children' not in cld4:
                            out_st_df += f'{n0}\t{n1}\t{n2}\t{n3}\t{n4}\n'
                        else:
                            for cld5 in cld4['children']:
                                n5 = cld5['name']
                                if 'children' not in cld5:
                                    out_st_df += f'{n0}\t{n1}\t{n2}\t{n3}\t{n4}\t{n5}\n'
                                else:
                                    for cld6 in cld5['children']:
                                        n6 = cld6['name']
                                        out_st_df += f'{n0}\t{n1}\t{n2}\t{n3}\t{n4}\t{n5}\t{n6}\n'
    return out_st_df
