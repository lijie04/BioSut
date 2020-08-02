
import os

def Version():
    ver1 = os.path.dirname(os.path.realpath(__file__)) + '/version'
    ver2 = os.path.dirname(os.path.realpath(__file__)) + '/../version'
    if os.path.isfile(ver1):return open(ver1).readline().strip()
    if os.path.isfile(ver2):return open(ver2).readline().strip()
    return 'Unknow version'
